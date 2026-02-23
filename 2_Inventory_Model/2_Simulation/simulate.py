# -*- coding: utf-8 -*-
"""Forward simulation / rolling evaluation of inventory policies.

Implements the empirical evaluation described in SS5.2 of:
  Mandl & Minner, "When Do Commodity Spot Price Regimes Matter for
  Inventory Managers?"

The simulation steps through actual monthly commodity prices. At each
period it calls the existing SDP solver to get the optimal first-period
ordering decision under each policy, executes the order at the realized
price, and tracks cumulative costs.

Usage
-----
  # Corn (non-Markov):
  python simulate.py --commodity corn [--n_months N]

  # Zinc (Markov):
  python simulate.py --commodity zinc [--n_months N]

Outputs a results table matching the format of paper Tables 2-3.

Runtime notes
-------------
- Corn (non-Markov): ~5-6 min for 104 months x 4 SDP policies x 3 ch scenarios
- Zinc (Markov):    ~30 min for 140 months x 4 SDP policies x 3 ch scenarios
  (Markov SDP is slower due to (last_price, current_price) Cartesian belief space)
- Use --n_months N to limit to the first N periods for quick testing
- Use --policies to run a subset of policies (e.g. --policies NAIVE PF MRS)
"""
import sys
import os
import argparse
import csv
import numpy as np

# ---------------------------------------------------------------------------
# Path setup: add Optimizer and Optimizer_Case directories
# ---------------------------------------------------------------------------
_SIM_DIR = os.path.dirname(os.path.abspath(__file__))
_INV_MODEL_DIR = os.path.join(_SIM_DIR, '..', '1_Model')
_OPTIMIZER_DIR = os.path.join(_INV_MODEL_DIR, 'Optimizer')
_CASE_DIR_BASE = os.path.join(_SIM_DIR, '..', '3_Case')

sys.path.insert(0, _OPTIMIZER_DIR)
sys.path.insert(0, os.path.join(_INV_MODEL_DIR, 'Optimizer_Case'))

from src_py.backward_recursion import backward_recursion
from src_py.backward_recursion_markov import backward_recursion_markov
from load_case import load_case_non_markov, load_case_markov

# Policy codes (for internal use)
_MRS    = 'MRS'
_CEC    = 'CEC'
_SRC_R1 = 'SRC-R1'
_SRC_R2 = 'SRC-R2'
_NAIVE  = 'NAIVE'
_PF     = 'PF'

# method_code values (must match load_case.py / update_pi.py)
_METHOD_EXACT = 0
_METHOD_CEC   = 1
_METHOD_R1    = 2
_METHOD_R2    = 3


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def load_historical_prices(csv_path):
    """Load monthly prices from a semicolon-delimited CSV (Date;Price).

    Returns (dates, prices) as lists.
    """
    dates, prices = [], []
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=';')
        for i, row in enumerate(reader):
            if i == 0 and row[0].strip().lower() in ('date', 'datum'):
                continue  # skip header
            if len(row) < 2 or not row[1].strip():
                continue
            dates.append(row[0].strip())
            prices.append(float(row[1].strip()))
    return dates, np.array(prices)


def snap_to_grid(price, price_grid):
    """Return the index of the nearest grid point to *price*."""
    return int(np.argmin(np.abs(price_grid - price)))


def bayesian_update_iid(pi, p_idx, k, P_pr_reg):
    """Simulation-level Bayesian belief update for i.i.d. regime switching.

    pi^j_{t+1} = sum_n( pi^n * k[n,j] * P_pr_reg[p_idx,n] )
                 / sum_n( pi^n * P_pr_reg[p_idx,n] )

    Parameters
    ----------
    pi : np.ndarray (m,)  current belief
    p_idx : int           index of the observed price on the grid
    k : np.ndarray (m,m)  transition matrix, rows sum to 1
    P_pr_reg : np.ndarray (N_prices, m)

    Returns
    -------
    pi_new : np.ndarray (m,)
    """
    m = len(pi)
    numerator = np.zeros(m)
    denominator = 0.0
    for n in range(m):
        denominator += pi[n] * P_pr_reg[p_idx, n]
    for j in range(m):
        for n in range(m):
            numerator[j] += pi[n] * k[n, j] * P_pr_reg[p_idx, n]
    pi_new = numerator / denominator
    return pi_new


def bayesian_update_markov(pi, p_idx, last_p_idx, k, R1, R2):
    """Simulation-level Bayesian belief update for Markov regime switching.

    pi^j_{t+1} = sum_n( pi^n * k[n,j] * R_n[last_p_idx, p_idx] )
                 / sum_n( pi^n * R_n[last_p_idx, p_idx] )

    Parameters
    ----------
    pi : np.ndarray (2,)
    p_idx, last_p_idx : int  current and previous price grid indices
    k : np.ndarray (2,2)
    R1, R2 : np.ndarray (N,N)  Gavirneni transition matrices

    Returns
    -------
    pi_new : np.ndarray (2,)
    """
    R = [R1, R2]
    m = len(pi)
    denominator = sum(pi[n] * R[n][last_p_idx, p_idx] for n in range(m))
    if denominator < 1e-15:
        # Unobserved transition - fall back to prior
        return pi.copy()
    pi_new = np.zeros(m)
    for j in range(m):
        for n in range(m):
            pi_new[j] += pi[n] * k[n, j] * R[n][last_p_idx, p_idx]
    pi_new /= denominator
    return pi_new


def solve_pf_full(future_prices, I_start, d, I_max, ch, cp):
    """Deterministic DP with known future prices (perfect foresight)."""
    n = len(future_prices)
    INF = 1e18

    # Forward DP storing optimal orders
    # V[t, I] = min cost from period t with inventory I
    V = np.full((n + 1, I_max + 1), 0.0)
    policy = np.zeros((n, I_max + 1), dtype=int)

    for t in range(n - 1, -1, -1):
        p_t = future_prices[t]
        for I in range(I_max + 1):
            best = INF
            for y in range(I_max - I + 1):
                net = I + y - d
                cost = p_t * y + (ch * net if net >= 0 else cp * (-net))
                I_next = max(net, 0)
                total = cost + V[t + 1, I_next]
                if total < best:
                    best = total
                    policy[t, I] = y
            V[t, I] = best

    return int(policy[0, I_start])


# ---------------------------------------------------------------------------
# Non-Markov simulation (Corn)
# ---------------------------------------------------------------------------

def simulate_non_markov(case_dir, historical_prices, ch_values, pi_0,
                        policies=None, n_months=None):
    """Rolling evaluation for the non-Markov (i.i.d.) case.

    Parameters
    ----------
    case_dir : str
        Path to the case data directory (Corn_Empirical).
    historical_prices : np.ndarray
        Monthly price series (actual realized prices).
    ch_values : list[float]
        Holding cost scenarios to evaluate.
    pi_0 : np.ndarray (m,)
        Initial regime belief.
    policies : list[str], optional
        Subset of policies to run. Default: all.
    n_months : int, optional
        Limit the number of months evaluated (for testing).

    Returns
    -------
    results : dict  {ch: {policy: {metric: value}}}
    """
    if policies is None:
        policies = [_MRS, _CEC, _SRC_R1, _SRC_R2, _NAIVE, _PF]

    # Load case parameters (ignore ch from config - overridden per scenario)
    (variables, _, sub_R1, sub_R2, sub_CEC,
     _, _, dist_code, _) = load_case_non_markov(case_dir)

    price_grid = variables['prices']
    k = variables['k']
    P_pr_reg = variables['P_pr_reg']
    d = variables['d_max']  # deterministic demand = 1 unit
    I_max = variables['I_max']
    t_max = variables['t_max']
    cp = variables['cp']

    T = len(historical_prices)
    if n_months is not None:
        T = min(T, n_months)

    results = {}

    for ch in ch_values:
        variables['ch'] = ch
        results[ch] = {}
        print(f"\n--- ch = {ch} ---")

        for policy in policies:
            print(f"  Running {policy}...", end=' ', flush=True)

            I = 0
            pi = pi_0.copy()
            total_purchase_cost = 0.0
            total_holding_cost = 0.0
            total_shortage_cost = 0.0
            total_units_purchased = 0
            inv_levels = []
            forward_buy_periods = 0

            for t in range(T):
                p_t = historical_prices[t]
                p_idx = snap_to_grid(p_t, price_grid)
                p_grid = price_grid[p_idx]

                # ---- Ordering decision ----
                if policy == _NAIVE:
                    y = max(d - I, 0)

                elif policy == _PF:
                    # Perfect foresight: look ahead t_max periods
                    horizon = min(t_max, T - t)
                    future_p = historical_prices[t:t + horizon]
                    y = solve_pf_full(future_p, I, d, I_max, ch, cp)

                else:
                    # SDP-based policies
                    if policy == _MRS:
                        method_code = _METHOD_EXACT
                    elif policy == _CEC:
                        method_code = _METHOD_CEC
                    elif policy == _SRC_R1:
                        method_code = _METHOD_R1
                    elif policy == _SRC_R2:
                        method_code = _METHOD_R2
                    else:
                        raise ValueError(f"Unknown policy: {policy}")

                    variables['pi'] = pi.copy()
                    variables['t_max'] = t_max

                    periods = backward_recursion(
                        variables, p_grid, I, dist_code,
                        method_code, 0, 0, 0  # dummy suboptimal decisions
                    )
                    y = int(periods[0]['period_min_order'][0])

                # ---- Realize costs ----
                net = I + y - d
                purchase_cost = p_t * y
                if net >= 0:
                    h_cost = ch * net
                    s_cost = 0.0
                else:
                    h_cost = 0.0
                    s_cost = cp * (-net)

                total_purchase_cost += purchase_cost
                total_holding_cost += h_cost
                total_shortage_cost += s_cost
                total_units_purchased += y
                inv_levels.append(max(net, 0))

                if y > max(d - I, 0):
                    forward_buy_periods += 1

                # ---- Update state ----
                I = max(net, 0)

                if policy == _MRS:
                    pi = bayesian_update_iid(pi, p_idx, k, P_pr_reg)

            total_cost = total_purchase_cost + total_holding_cost + total_shortage_cost
            avg_price_paid = (total_purchase_cost / total_units_purchased
                              if total_units_purchased > 0 else 0.0)
            avg_inventory = np.mean(inv_levels)

            results[ch][policy] = {
                'total_cost': total_cost,
                'total_purchase': total_purchase_cost,
                'total_holding': total_holding_cost,
                'total_shortage': total_shortage_cost,
                'avg_price_paid': avg_price_paid,
                'avg_inventory': avg_inventory,
                'forward_buy_periods': forward_buy_periods,
                'n_periods': T,
            }
            print(f"total_cost={total_cost:.2f}")

    return results


# ---------------------------------------------------------------------------
# Markov simulation (Zinc)
# ---------------------------------------------------------------------------

def simulate_markov(case_dir, historical_prices, ch_values, pi_0,
                    policies=None, n_months=None):
    """Rolling evaluation for the Markov case.

    Parameters
    ----------
    case_dir : str
        Path to the case data directory (Zinc).
    historical_prices : np.ndarray
    ch_values : list[float]
    pi_0 : np.ndarray (2,)
    policies : list[str], optional
    n_months : int, optional

    Returns
    -------
    results : dict  {ch: {policy: {metric: value}}}
    """
    if policies is None:
        policies = [_MRS, _CEC, _SRC_R1, _SRC_R2, _NAIVE, _PF]

    (variables, _, sub_R1, sub_R2, sub_CEC,
     _, _, dist_code, _) = load_case_markov(case_dir)

    price_grid = variables['prices']
    k = variables['k']
    R1 = variables['R1']
    R2 = variables['R2']
    d = variables['d_max']
    I_max = variables['I_max']
    t_max = variables['t_max']
    cp = variables['cp']

    T = len(historical_prices)
    if n_months is not None:
        T = min(T, n_months)

    price_dist = ['R1', 'R2']
    results = {}

    for ch in ch_values:
        variables['ch'] = ch
        results[ch] = {}
        print(f"\n--- ch = {ch} ---")

        for policy in policies:
            print(f"  Running {policy}...", end=' ', flush=True)

            I = 0
            pi = pi_0.copy()
            # Markov case needs last_price: initialise to first observation
            last_p_idx = snap_to_grid(historical_prices[0], price_grid)

            total_purchase_cost = 0.0
            total_holding_cost = 0.0
            total_shortage_cost = 0.0
            total_units_purchased = 0
            inv_levels = []
            forward_buy_periods = 0

            for t in range(T):
                p_t = historical_prices[t]
                p_idx = snap_to_grid(p_t, price_grid)
                p_grid = price_grid[p_idx]
                last_p_grid = price_grid[last_p_idx]

                # ---- Ordering decision ----
                if policy == _NAIVE:
                    y = max(d - I, 0)

                elif policy == _PF:
                    horizon = min(t_max, T - t)
                    future_p = historical_prices[t:t + horizon]
                    y = solve_pf_full(future_p, I, d, I_max, ch, cp)

                else:
                    if policy == _MRS:
                        method_code = _METHOD_EXACT
                    elif policy == _CEC:
                        method_code = _METHOD_CEC
                    elif policy == _SRC_R1:
                        method_code = _METHOD_R1
                    elif policy == _SRC_R2:
                        method_code = _METHOD_R2
                    else:
                        raise ValueError(f"Unknown policy: {policy}")

                    variables['pi'] = pi.copy()
                    variables['t_max'] = t_max

                    periods = backward_recursion_markov(
                        variables, p_grid, last_p_grid, I, dist_code,
                        price_dist, method_code, 0, 0, 0
                    )
                    if periods is None:
                        # No reachable state - fall back to NAIVE
                        y = max(d - I, 0)
                    else:
                        y = int(periods[0]['period_min_order'][0])

                # ---- Realize costs ----
                net = I + y - d
                purchase_cost = p_t * y
                if net >= 0:
                    h_cost = ch * net
                    s_cost = 0.0
                else:
                    h_cost = 0.0
                    s_cost = cp * (-net)

                total_purchase_cost += purchase_cost
                total_holding_cost += h_cost
                total_shortage_cost += s_cost
                total_units_purchased += y
                inv_levels.append(max(net, 0))

                if y > max(d - I, 0):
                    forward_buy_periods += 1

                # ---- Update state ----
                prev_p_idx = p_idx
                I = max(net, 0)

                if policy == _MRS:
                    pi = bayesian_update_markov(pi, p_idx, last_p_idx, k, R1, R2)

                last_p_idx = prev_p_idx

            total_cost = total_purchase_cost + total_holding_cost + total_shortage_cost
            avg_price_paid = (total_purchase_cost / total_units_purchased
                              if total_units_purchased > 0 else 0.0)
            avg_inventory = np.mean(inv_levels)

            results[ch][policy] = {
                'total_cost': total_cost,
                'total_purchase': total_purchase_cost,
                'total_holding': total_holding_cost,
                'total_shortage': total_shortage_cost,
                'avg_price_paid': avg_price_paid,
                'avg_inventory': avg_inventory,
                'forward_buy_periods': forward_buy_periods,
                'n_periods': T,
            }
            print(f"total_cost={total_cost:.2f}")

    return results


# ---------------------------------------------------------------------------
# Results display
# ---------------------------------------------------------------------------

def print_results_table(results, commodity):
    """Print results in the style of paper Tables 2-3."""
    print(f"\n{'='*70}")
    print(f"  {commodity.upper()} - Rolling Evaluation Results")
    print(f"{'='*70}")

    for ch, policy_results in results.items():
        print(f"\n  Holding cost ch = {ch}")
        print(f"  {'Policy':<10} {'Total Cost':>14} {'% > PF':>8} "
              f"{'Avg Price':>11} {'Avg Inv':>9} {'Fwd Buy':>8}")
        print(f"  {'-'*65}")

        pf_cost = policy_results.get(_PF, {}).get('total_cost', None)

        for policy, m in policy_results.items():
            tc = m['total_cost']
            pct_above_pf = ((tc / pf_cost - 1) * 100
                            if pf_cost and pf_cost > 0 else float('nan'))
            print(f"  {policy:<10} {tc:>14.2f} {pct_above_pf:>7.1f}% "
                  f"{m['avg_price_paid']:>11.2f} {m['avg_inventory']:>9.3f} "
                  f"{m['forward_buy_periods']:>8d}")

    print(f"\n{'='*70}\n")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Rolling evaluation of inventory policies on historical price data.')
    parser.add_argument('--commodity', choices=['corn', 'zinc'], default='corn',
                        help='Which commodity to evaluate (default: corn)')
    parser.add_argument('--n_months', type=int, default=None,
                        help='Limit simulation to first N months (for testing)')
    parser.add_argument('--policies', nargs='+',
                        choices=[_MRS, _CEC, _SRC_R1, _SRC_R2, _NAIVE, _PF],
                        default=None,
                        help='Policies to run (default: all)')
    args = parser.parse_args()

    if args.commodity == 'corn':
        case_dir = os.path.join(_CASE_DIR_BASE, 'Corn_Empirical')
        csv_path = os.path.join(_SIM_DIR, '..', '..', '1_Price_Model',
                                '3_ARIMA', 'csv_Input', 'Corn_m_red.csv')
        # ch in cents/bu/month: paper values {0.01, 0.10, 0.40} $/bu = {1, 10, 40} cents/bu
        ch_values = [1.0, 10.0, 40.0]
        pi_0 = np.array([0.375, 0.625])

        _, hist_prices = load_historical_prices(csv_path)
        print(f"Corn: {len(hist_prices)} monthly observations loaded.")

        results = simulate_non_markov(
            case_dir, hist_prices, ch_values, pi_0,
            policies=args.policies, n_months=args.n_months
        )
        print_results_table(results, 'Corn')

    else:  # zinc
        case_dir = os.path.join(_CASE_DIR_BASE, 'Zinc')
        csv_path = os.path.join(_SIM_DIR, '..', '..', '1_Price_Model',
                                '3_ARIMA', 'csv_Input', 'Zinc_m.csv')
        ch_values = [17.5, 87.5, 350.0]  # examples; paper uses multiple scenarios
        pi_0 = np.array([0.5, 0.5])

        _, hist_prices = load_historical_prices(csv_path)
        print(f"Zinc: {len(hist_prices)} monthly observations loaded.")

        results = simulate_markov(
            case_dir, hist_prices, ch_values, pi_0,
            policies=args.policies, n_months=args.n_months
        )
        print_results_table(results, 'Zinc')


if __name__ == '__main__':
    main()
