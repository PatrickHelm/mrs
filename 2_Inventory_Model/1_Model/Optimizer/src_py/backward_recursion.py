import numpy as np
from numba import njit

from .single_period_cost import single_period_cost
from .backward_cost_function import backward_cost_function
from .get_demand_dist import get_demand_dist
from .update_pi import update_pi, OLFC
from .next_price import next_price
from .needed_pi import needed_pi
from .get_next_prices import get_next_prices
from .reorder import aux_matrix
from .get_mins import get_mins, precompute_mins_indices


def combvec(*arrays):
    """Create all combinations (equivalent to MATLAB combvec transposed).

    First argument varies fastest, matching MATLAB's combvec behavior.
    """
    grids = np.meshgrid(*arrays, indexing='ij')
    return np.column_stack([g.ravel(order='F') for g in grids])


@njit(cache=True)
def _compute_period_last(all_comb, I_max, cp, ch, d_max, demand_dist, dist_code):
    """JIT-compiled s×y loop for the last period (terminal).

    Only single_period_cost is needed (no future costs).
    """
    n_states = all_comb.shape[0]
    period_matrix = np.zeros((n_states, I_max + 1))
    period_min_cost = np.zeros(n_states)
    period_min_order = np.zeros(n_states, dtype=np.int64)

    for s in range(n_states):
        p = all_comb[s, 0]
        I = int(all_comb[s, 1])
        for y in range(I_max + 1):
            period_matrix[s, y] = single_period_cost(
                y, cp, ch, p, d_max, I, demand_dist, dist_code
            )
        period_min_cost[s] = np.min(period_matrix[s])
        period_min_order[s] = np.argmin(period_matrix[s])

    return period_matrix, period_min_cost, period_min_order


@njit(cache=True)
def _compute_period(all_comb, I_max, cp, ch, d_max, demand_dist, dist_code,
                    next_prices_2d, mins_flat, mins_indices):
    """JIT-compiled s×y loop for non-terminal periods.

    Uses precomputed mins_flat/mins_indices and next_prices_2d.
    """
    n_states = all_comb.shape[0]
    period_matrix = np.zeros((n_states, I_max + 1))
    period_min_cost = np.zeros(n_states)
    period_min_order = np.zeros(n_states, dtype=np.int64)

    for s in range(n_states):
        p = all_comb[s, 0]
        I = int(all_comb[s, 1])
        next_price_vec = next_prices_2d[s]
        mins = mins_flat[mins_indices[s, 0]:mins_indices[s, 1]]

        for y in range(I_max + 1):
            period_matrix[s, y] = backward_cost_function(
                y, cp, ch, p, d_max, I_max, I,
                next_price_vec, mins, demand_dist, dist_code
            )
        period_min_cost[s] = np.min(period_matrix[s])
        period_min_order[s] = np.argmin(period_matrix[s])

    return period_matrix, period_min_cost, period_min_order


def backward_recursion(variables, observed_price, start_inventory, dist_code,
                       method_code, suboptimal_decision_CEC,
                       suboptimal_decision_R1, suboptimal_decision_R2):
    """Backward dynamic programming recursion (non-Markov).

    Parameters
    ----------
    variables : dict
        Must contain keys: ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices,
        and distribution parameters (normal_mu, normal_sigma, poisson_lambda, nbin_r, nbin_p).
    observed_price : float
    start_inventory : int
    dist_code : int
        0=uniform, 1=normal, 2=poisson, 3=nbin, 4=manual.
    method_code : int
        0=exact, 1=CEC, 2=R1, 3=R2, 4=OLFC.
    suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2 : int

    Returns
    -------
    periods : list of dict
        One dict per period with keys: PI, PI_needed, period_matrix,
        period_min_cost, period_min_order, suboptimal_cost_CEC/R1/R2, next_price.
    """
    # Unpack variables
    ch = variables['ch']
    cp = variables['cp']
    d_max = variables['d_max']
    I_max = variables['I_max']
    k = variables['k']
    m = variables['m']
    pi = variables['pi']
    P_pr_reg = variables['P_pr_reg']
    t_max = variables['t_max']
    prices = variables['prices']

    # Get demand distribution
    demand_dist = get_demand_dist(
        dist_code, d_max,
        normal_mu=variables.get('normal_mu', 15.0),
        normal_sigma=variables.get('normal_sigma', 3.0),
        poisson_lambda=variables.get('poisson_lambda', 1.0),
        nbin_r=variables.get('nbin_r', 22.5),
        nbin_p=variables.get('nbin_p', 0.6),
        manual_dist=variables.get('manual_dist', None),
    )

    # Ensure pi is 2D: (1, m)
    pi_2d = np.atleast_2d(pi).astype(np.float64)

    # --- Forward pass: compute beliefs and next-price probabilities ---
    PI = [None] * t_max
    nP = [None] * t_max

    for t in range(t_max):
        if t == 0:
            PI[t] = pi_2d.copy()
        else:
            if method_code == OLFC:
                PI[t] = update_pi(k, m, pi_2d, P_pr_reg, t + 1, method_code)
            else:
                PI[t] = update_pi(k, m, PI[t - 1], P_pr_reg, t + 1, method_code)
        nP[t] = next_price(PI[t], P_pr_reg)

    # Filter to reachable beliefs
    PI_needed = needed_pi(PI, observed_price, prices)

    inv = np.arange(I_max + 1, dtype=np.float64)

    # --- Backward pass ---
    periods = [None] * t_max

    for t in range(t_max - 1, -1, -1):
        t_matlab = t + 1  # 1-based period number

        if t == t_max - 1:
            # Last period: only single-period costs
            all_comb = combvec(prices, inv, PI_needed[t])
            period_matrix, period_min_cost, period_min_order = _compute_period_last(
                all_comb, I_max, cp, ch, d_max, demand_dist, dist_code
            )
        else:
            # Non-last period
            if t != 0:
                all_comb = combvec(prices, inv, PI_needed[t])
            else:
                all_comb = combvec(
                    np.array([observed_price]),
                    np.array([float(start_inventory)]),
                    PI_needed[t]
                )

            n_states = all_comb.shape[0]

            # Get next-period price probabilities
            next_prices_raw = get_next_prices(
                nP[t + 1], all_comb[:, 0], observed_price, prices, t_matlab
            )

            if t != 0:
                next_prices_2d = aux_matrix(
                    next_prices_raw, len(prices), I_max, t_matlab
                )
            else:
                next_prices_2d = next_prices_raw

            # Ensure next_prices_2d matches n_states
            if next_prices_2d.shape[0] < n_states:
                # For t==0, only one state; replicate if needed
                next_prices_2d = np.tile(next_prices_2d, (n_states, 1))

            # Precompute mins for all states
            future_min_cost = periods[t + 1]['period_min_cost']
            mins_flat, mins_indices = precompute_mins_indices(
                PI_needed, prices, t, future_min_cost, inv, all_comb
            )

            period_matrix, period_min_cost, period_min_order = _compute_period(
                all_comb, I_max, cp, ch, d_max, demand_dist, dist_code,
                next_prices_2d, mins_flat, mins_indices
            )

        # Store results
        n_states = all_comb.shape[0]
        s_last = n_states - 1  # MATLAB s retains last loop value

        periods[t] = {
            'PI': PI[t],
            'PI_needed': PI_needed[t],
            'period_matrix': period_matrix,
            'suboptimal_cost_CEC': period_matrix[s_last, suboptimal_decision_CEC],
            'suboptimal_cost_R1': period_matrix[s_last, suboptimal_decision_R1],
            'suboptimal_cost_R2': period_matrix[s_last, suboptimal_decision_R2],
            'period_min_cost': period_min_cost,
            'period_min_order': period_min_order,
            'next_price': nP[t],
        }

    return periods
