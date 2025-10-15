import numpy as np

def backward_cost_function(y, cp, ch, p, d_max, I_max, I, next_price, mins, demand_dist, dist):
    d_clean = d_max + 1  # 0 can also be an observed demand
    cost = 0
    period_cost = single_period_cost(y, cp, ch, p, d_max, I, demand_dist, dist)  # Costs if t would be the last period

    # Probability that Inventory gets 0
    probIgets0 = np.sum(demand_dist[I + y : d_clean])

    # Probability that inventory becomes > 0 -> [d=m ... d=2 d=1 d=0]
    probIgreater0 = []
    if I + y <= d_max:
        for d in range(I + y - 1, -1, -1):
            probIgreater0.append(demand_dist[d])
    else:
        probIgreater0 = demand_dist.copy()

    probIgreater0 = np.array(probIgreater0)
    prob_price_mat = np.outer(probIgreater0, next_price)  # Matrix of probabilities * next price probabilities

    block = len(next_price)  # Number of different prices
    ende = block * (I + y + 1)  # Last value of all minimum costs from t+1 needed for calculation

    if I + y == 0:  # Inventory in t+1 is always 0 -> probIgets0=1;
        for i in range(ende):
            cost += probIgets0 * next_price[i % block] * mins[i]
    elif 0 < I + y <= I_max and I + y - d_max <= 0:  # I+y cannot be > I_max
        for i in range(block):
            cost += probIgets0 * next_price[i] * mins[i]
        j = block
        for P in range(len(probIgreater0)):  # Number of possible scenarios where inventory > 0 in t+1
            for np_idx in range(len(next_price)):
                cost += prob_price_mat[P, np_idx] * mins[j]
                j += 1
    else:  # I+y > I_max
        j = (I + y - d_max) * block
        for d in range(d_clean, 0, -1):
            if I + y - (d - 1) <= I_max:
                for np_idx in range(len(next_price)):
                    cost += prob_price_mat[d - 1, np_idx] * mins[j]
                    j += 1
            else:
                cost += period_cost
                return cost

    cost += period_cost
    return cost