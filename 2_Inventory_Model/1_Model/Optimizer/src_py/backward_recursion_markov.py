import numpy as np
from numba import njit

from .single_period_cost import single_period_cost
from .backward_cost_function import backward_cost_function
from .get_demand_dist import get_demand_dist
from .get_markov_model import get_markov_model
from .update_pi_markov import update_pi_markov
from .next_price_markov import next_price_markov
from .needed_pi_markov import needed_pi_markov
from .get_next_prices_markov import get_next_prices_markov
from .reorder_markov import reorder_markov
from .get_mins import get_mins, precompute_mins_indices


def combvec(*arrays):
    """Create all combinations (equivalent to MATLAB combvec transposed).

    First argument varies fastest, matching MATLAB's combvec behavior.
    """
    grids = np.meshgrid(*arrays, indexing='ij')
    return np.column_stack([g.ravel(order='F') for g in grids])


@njit(cache=True)
def _compute_period_last(all_comb, I_max, cp, ch, d_max, demand_dist, dist_code):
    """JIT-compiled s×y loop for the last period."""
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
    """JIT-compiled s×y loop for non-terminal periods."""
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


def backward_recursion_markov(variables, observed_price, last_price,
                              start_inventory, dist_code, price_dist,
                              method_code, suboptimal_decision_CEC,
                              suboptimal_decision_R1, suboptimal_decision_R2):
    """Backward dynamic programming recursion (Markov version).

    Parameters
    ----------
    variables : dict
        Must contain keys: ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices,
        and distribution parameters.
    observed_price : float
    last_price : float
    start_inventory : int
    dist_code : int
        0=uniform, 1=normal, 2=poisson, 3=nbin, 4=manual.
    price_dist : list of str
        e.g. ['MR', 'MO']
    method_code : int
        0=exact, 1=CEC, 2=R1, 3=R2.
    suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2 : int

    Returns
    -------
    periods : list of dict, or None if no reachable state.
    """
    # Unpack variables
    ch = variables['ch']
    cp = variables['cp']
    d_max = variables['d_max']
    I_max = variables['I_max']
    k = variables['k']
    m = variables['m']
    pi = variables['pi']
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

    # Build Markov model (use pre-built if provided in variables)
    if 'R1' in variables and 'R2' in variables:
        R1, R2 = variables['R1'], variables['R2']
    else:
        R1, R2 = get_markov_model(price_dist)

    # --- Forward pass: compute Markov beliefs and next-price probabilities ---
    # MATLAB: for t=1:t_max+1 (builds t_max+1 entries, then slices)
    PI_markov_raw = [None] * (t_max + 1)
    nP_markov_raw = [None] * (t_max + 1)

    pi_arr = np.array(pi, dtype=np.float64)

    for t in range(t_max + 1):
        t_matlab = t + 1  # 1-based
        if t == 0:
            PI_markov_raw[t] = update_pi_markov(R1, R2, prices, t_matlab, pi_arr, k, method_code)
            nP_markov_raw[t] = next_price_markov(R1, R2, prices, t_matlab, pi_arr)
        else:
            PI_markov_raw[t] = update_pi_markov(R1, R2, prices, t_matlab, PI_markov_raw, k, method_code)
            nP_markov_raw[t] = next_price_markov(R1, R2, prices, t_matlab, PI_markov_raw)

    # MATLAB: PI_markov = PI_markov(2:end); nP_markov = nP_markov(1:end-1)
    PI_markov = PI_markov_raw[1:]      # periods 1..t_max (0-indexed: 0..t_max-1)
    nP_markov = nP_markov_raw[:t_max]  # periods 1..t_max

    # Filter reachable beliefs
    PI_needed_markov = needed_pi_markov(PI_markov, observed_price, last_price, prices)

    # Check for no reachable state
    if (PI_needed_markov[0].shape[0] == 1
            and np.allclose(PI_needed_markov[0], np.array([[0.0, 0.0]]))):
        print('No reachable state')
        return None

    # Extract column 0 (R1 beliefs) as PI_needed for get_mins compatibility
    PI_needed = [None] * t_max
    for t in range(t_max):
        PI_needed[t] = PI_needed_markov[t][:, 0]

    inv = np.arange(I_max + 1, dtype=np.float64)

    # --- Backward pass ---
    periods = [None] * t_max

    for t in range(t_max - 1, -1, -1):
        t_matlab = t + 1

        if t == t_max - 1:
            # Last period
            all_comb = combvec(prices, inv, PI_needed[t])
            period_matrix, period_min_cost, period_min_order = _compute_period_last(
                all_comb, I_max, cp, ch, d_max, demand_dist, dist_code
            )
        else:
            if t != 0:
                all_comb = combvec(prices, inv, PI_needed[t])
            else:
                all_comb = combvec(
                    np.array([observed_price]),
                    np.array([float(start_inventory)]),
                    PI_needed[t]
                )

            n_states = all_comb.shape[0]

            # Get next-period price probabilities (Markov)
            next_prices_raw = get_next_prices_markov(
                nP_markov[t + 1], all_comb[:, 0],
                observed_price, last_price, prices, t_matlab
            )

            if t != 0:
                next_prices_2d = reorder_markov(
                    next_prices_raw, len(prices), I_max, t_matlab
                )
            else:
                next_prices_2d = next_prices_raw

            if next_prices_2d.shape[0] < n_states:
                next_prices_2d = np.tile(next_prices_2d, (n_states, 1))

            # Precompute mins
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
        s_last = n_states - 1

        periods[t] = {
            'PI': PI_markov[t],
            'PI_needed': PI_needed[t],
            'period_matrix': period_matrix,
            'suboptimal_cost_CEC': period_matrix[s_last, suboptimal_decision_CEC],
            'suboptimal_cost_R1': period_matrix[s_last, suboptimal_decision_R1],
            'suboptimal_cost_R2': period_matrix[s_last, suboptimal_decision_R2],
            'period_min_cost': period_min_cost,
            'period_min_order': period_min_order,
            'nextPrice': nP_markov[t],
        }

    return periods
