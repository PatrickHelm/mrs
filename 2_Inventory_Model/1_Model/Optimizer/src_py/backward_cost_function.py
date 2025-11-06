import numpy as np
from .single_period_cost import single_period_cost

def backward_cost_function(
    procurement_decision: int,
    penalty_cost: float,
    holding_cost: float,
    price: float,
    d_max: int,
    I_max: int,
    inventory: int,
    next_price: np.ndarray,   # 1D array, length = number of price states
    min_value_functions: np.ndarray,         # 1D array of future value function (min costs)
    demand_dist: np.ndarray,  # 1D array, P(D = 0..d_max)
) -> float:
    I_y = inventory + procurement_decision
    block = next_price.shape[0]

    # Cost if current period were the last one
    period_cost = single_period_cost(
        procurement_decision,
        penalty_cost,
        holding_cost,
        price,
        d_max,
        inventory,
        demand_dist,
    )

    # Probability that inventory in t+1 is zero: sum_{d >= I_y} P(D = d)
    prob_inv_zero = demand_dist[I_y: d_max + 1].sum()

    cost = 0.0

    # --- Case 1: I + y == 0 ---------------------------------------------
    if I_y == 0:
        # ende = block * (I_y + 1) = block
        cost += prob_inv_zero * np.dot(next_price, min_value_functions[:block])

    # --- Case 2: 0 < I + y <= I_max and I + y <= d_max -------------------
    elif (I_y > 0) and (I_y <= I_max) and (I_y - d_max <= 0):
        # Probabilities that inventory becomes > 0:
        # [P(D = I_y - 1), ..., P(D = 0)]
        prob_inv_gt_zero = demand_dist[:I_y][::-1]  # length I_y

        # Matrix of joint probabilities (positive inventory, next price)
        prob_price_mat = prob_inv_gt_zero[:, None] * next_price[None, :]

        # First block: inventory becomes 0 in t+1
        cost += prob_inv_zero * np.dot(next_price, min_value_functions[:block])

        # Remaining blocks: inventory > 0 in t+1
        n_pos = prob_inv_gt_zero.shape[0]
        mins_block = min_value_functions[block : block + n_pos * block].reshape(n_pos, block)
        cost += np.sum(prob_price_mat * mins_block)

    # --- Case 3: else branch (like MATLAB 'else % I+y>I_max') ------------
    else:
        # Here MATLAB uses probIgreater0 = demand_dist (full vector)
        prob_inv_gt_zero_full = demand_dist
        prob_price_mat_full = prob_inv_gt_zero_full[:, None] * next_price[None, :]

        # j = (I + y - d_max) * block + 1   (MATLAB, 1-based)
        j = (I_y - d_max) * block          # Python, 0-based

        # d loops from d_max down to 0 (actual demand)
        for d in range(d_max, -1, -1):
            if I_y - d <= I_max:
                for np_idx in range(block):
                    cost += prob_price_mat_full[d, np_idx] * min_value_functions[j]
                    j += 1
            else:
                cost += period_cost
                return cost + period_cost

    # Add current-period cost
    return cost + period_cost