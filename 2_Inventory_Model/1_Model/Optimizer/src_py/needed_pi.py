import numpy as np
from itertools import product

def needed_pi(PI, observed_price, prices):
    t = len(PI)
    res = [None] * t

    actual_price = np.array(prices) == observed_price

    for i in range(t):
        if i == 0:  # Reachable belief in period t=1
            res[0] = PI[0][0]
        elif i == 1:  # Reachable belief in period t=2
            res[1] = PI[1][actual_price, 0]
        else:
            counter = prices
            for j in range(i - 1):
                var = list(product(counter, prices))
                counter = var
                var = np.array(var)[:, ::-1]
                ind = var[:, 0] == observed_price
                res[i] = PI[i][ind, 0]

    return res
