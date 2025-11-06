import numpy as np

def reorder(next_prices, prices, I_max, t):
    I_max = I_max + 1  # labeling in Matlab
    num_prices = len(prices)  # number of different prices
    num_inv_price_combinations = I_max * num_prices  # number of inventory-price combinations

    res = np.zeros((num_prices**(t-1)*I_max, num_prices))
    start = 0
    a = 0
    b = num_inv_price_combinations
    while a < len(res):
        ende = start + (num_prices - 1)  # start:ende = rows that will be copied
        res[a:b, :] = np.tile(next_prices[start:ende+1, :], (I_max, 1))  # copies the rows times the number of possible inventory states
        start = ende + 1
        a = b
        b = b + num_inv_price_combinations
    return res