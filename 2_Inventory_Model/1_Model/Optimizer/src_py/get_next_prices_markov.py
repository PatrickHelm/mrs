import numpy as np
from itertools import product

def get_next_prices_markov(next_price, p, observed_price, last_price, prices, t):
    if t >= 1:
        var = prices
        for _ in range(t):
            var = np.array(list(product(var, prices)))
        var = np.fliplr(var)  # Creation of all combinations of P1, P2, ..., PN
        pInd1 = (var[:, 0] == last_price) & (var[:, 1] == observed_price)
        res = next_price[pInd1, :]
    else:  # t=1
        for i, price in enumerate(prices):
            if last_price == price:
                res = next_price[i, :]
                break
    return res