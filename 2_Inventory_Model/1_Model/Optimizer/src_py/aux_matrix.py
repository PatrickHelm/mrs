import numpy as np

def aux_matrix(next_prices, prices, I_max, t):
    I_max = I_max + 1  # labeling in Matlab
    fac1 = len(prices)  # number of different prices
    fac2 = I_max * fac1  # number of inventory-price combinations

    res = np.zeros((fac1**(t-1)*I_max, fac1))
    start = 0
    a = 0
    b = fac2
    while a < len(res):
        ende = start + (fac1 - 1)  # start:ende = rows that will be copied
        res[a:b, :] = np.tile(next_prices[start:ende+1, :], (I_max, 1))  # copies the rows times the number of possible inventory states
        start = ende + 1
        a = b
        b = b + fac2
    return res