import numpy as np
from numba import njit

# @njit
def single_period_cost(procurement_decision: int, 
                       penalty_cost: float, 
                       holding_cost: float, 
                       price: float, 
                       d_max: int, 
                       inventory: int, 
                       demand_dist: np.ndarray) -> float:
    # start with procurement cost
    costs = price * procurement_decision
    for d in range(d_max + 1):
        dProb = demand_dist[d]
        # add holding cost
        costs += holding_cost * dProb * max(0, inventory + procurement_decision - d)
        # add penalty cost
        costs += penalty_cost * dProb * max(0, d - inventory - procurement_decision)
    return costs