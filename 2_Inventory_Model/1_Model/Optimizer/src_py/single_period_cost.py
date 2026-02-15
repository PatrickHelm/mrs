import numpy as np
from numba import njit

DIST_UNIFORM = 0


@njit(cache=True)
def single_period_cost(y, cp, ch, p, d_max, I, demand_dist, dist_code):
    """Calculate the cost for a single (last) period given an order quantity.

    Parameters
    ----------
    y : int
        Order quantity.
    cp : float
        Shortage penalty cost per unit.
    ch : float
        Holding cost per unit.
    p : float
        Current price.
    d_max : int
        Maximum demand.
    I : int
        Current inventory level.
    demand_dist : np.ndarray
        Demand probability distribution (length d_max+1).
    dist_code : int
        0=uniform, otherwise general.

    Returns
    -------
    costs : float
    """
    costs = p * y

    if dist_code == DIST_UNIFORM:
        inv_prob = 1.0 / (d_max + 1)
        for d in range(d_max + 1):
            excess = I + y - d
            shortage = d - I - y
            costs += ch * inv_prob * max(excess, 0.0) + cp * inv_prob * max(shortage, 0.0)
    else:
        for d in range(d_max + 1):
            d_prob = demand_dist[d]
            excess = I + y - d
            shortage = d - I - y
            costs += ch * d_prob * max(excess, 0.0) + cp * d_prob * max(shortage, 0.0)

    return costs
