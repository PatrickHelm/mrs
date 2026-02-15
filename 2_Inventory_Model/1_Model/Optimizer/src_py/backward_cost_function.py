import numpy as np
from numba import njit
from .single_period_cost import single_period_cost

DIST_UNIFORM = 0


@njit(cache=True)
def backward_cost_function(y, cp, ch, p, d_max, I_max, I, next_price_vec,
                           mins, demand_dist, dist_code):
    """Calculate total expected cost for ordering y in a non-terminal period.

    Combines immediate period cost with expected future costs over all
    possible demand realizations and next-period prices.

    Parameters
    ----------
    y : int
        Order quantity.
    cp : float
        Shortage penalty cost.
    ch : float
        Holding cost.
    p : float
        Current price.
    d_max : int
        Maximum demand.
    I_max : int
        Maximum inventory.
    I : int
        Current inventory level.
    next_price_vec : np.ndarray, shape (n_prices,)
        Probabilities of each next-period price.
    mins : np.ndarray
        Minimum future costs from period t+1 (precomputed slice).
    demand_dist : np.ndarray, shape (d_max+1,)
        Demand probability distribution.
    dist_code : int
        0=uniform, else general.

    Returns
    -------
    cost : float
    """
    d_clean = d_max + 1  # number of demand levels (0 is also possible)
    cost = 0.0
    period_cost = single_period_cost(y, cp, ch, p, d_max, I, demand_dist, dist_code)

    # Probability that inventory becomes 0 (demand >= I+y)
    # MATLAB: sum(demand_dist(I+y+1:d_clean)) → Python: sum(demand_dist[I+y:])
    prob_I_gets_0 = 0.0
    for idx in range(I + y, d_clean):
        prob_I_gets_0 += demand_dist[idx]

    # Probability that inventory stays > 0: for d = I+y-1 down to 0
    # Result is ordered [d=I+y-1, d=I+y-2, ..., d=1, d=0]
    block = len(next_price_vec)

    if I + y <= d_max:
        n_prob = I + y  # number of demand values that leave inventory > 0
        prob_I_greater_0 = np.empty(n_prob)
        idx = 0
        for d in range(I + y - 1, -1, -1):
            prob_I_greater_0[idx] = demand_dist[d]
            idx += 1
    else:
        # All demands leave some inventory
        prob_I_greater_0 = demand_dist.copy()

    # prob_price_mat = probIgreater0' * next_price (outer product)
    n_prob_ig = len(prob_I_greater_0)
    prob_price_mat = np.empty((n_prob_ig, block))
    for pi in range(n_prob_ig):
        for npi in range(block):
            prob_price_mat[pi, npi] = prob_I_greater_0[pi] * next_price_vec[npi]

    # --- Three branches based on I+y ---

    if I + y == 0:
        # Inventory in t+1 is always 0 → probIgets0 = 1
        # MATLAB: ende = block*(I+y+1) = block*1 = block
        ende = block
        for i in range(ende):
            cost += prob_I_gets_0 * next_price_vec[i] * mins[i]

    elif I + y > 0 and I + y <= I_max and I + y - d_max <= 0:
        # Normal case
        # First: cost contribution from inventory becoming 0
        for i in range(block):
            cost += prob_I_gets_0 * next_price_vec[i] * mins[i]

        # Then: cost contributions from inventory staying > 0
        j = block  # continues from where the first loop left off
        for P in range(n_prob_ig):
            for npi in range(block):
                cost += prob_price_mat[P, npi] * mins[j]
                j += 1

    else:
        # I+y > I_max or I+y - d_max > 0: clipped case
        # MATLAB: j = (I+y-d_max)*block + 1 → Python: j = (I+y - d_max)*block
        j = (I + y - d_max) * block

        # MATLAB: for d = d_clean:-1:1 (d is 1-based index into demand levels)
        for d_idx in range(d_clean - 1, -1, -1):
            # MATLAB: d goes from d_clean down to 1
            # MATLAB: I+y-(d-1) → in Python: I+y - d_idx (since d_idx = d_matlab - 1)
            if I + y - d_idx <= I_max:
                for npi in range(block):
                    cost += prob_price_mat[d_idx, npi] * mins[j]
                    j += 1
            else:
                cost += period_cost
                return cost

    cost += period_cost
    return cost
