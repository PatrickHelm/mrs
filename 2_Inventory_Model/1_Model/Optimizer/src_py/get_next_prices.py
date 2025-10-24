import numpy as np
from itertools import product

def get_next_prices(next_price, p, observed_price, prices, t):
    if t >= 2:
        var = prices
        for _ in range(t - 1):
            var = np.array(list(product(var, prices)))
        var = np.fliplr(var)  # Creation of all combinations of P1, P2, ..., PN
        pInd1 = var[:, 0] == observed_price
        res = next_price[pInd1, :]
    else:  # t=1
        for i in range(len(prices)):
            if p == prices[i]:
                res = next_price[i, :]
                break
    return res