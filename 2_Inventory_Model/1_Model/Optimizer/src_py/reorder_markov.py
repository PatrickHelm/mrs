import numpy as np

def reorder_markov(next_prices, b_t, I_max, t):
    """
    Reorder creates a large matrix with copies of 'next_prices' to fit the
    observed states in periods t.
    """
    I_max += 1  # labelimg adjustment
    fac1 = len(b_t)                # Number of different prices
    fac2 = I_max * fac1           # Number of inventory states for each price

    num_rows = (fac1 ** t) * I_max
    res = np.zeros((num_rows, fac1))
    
    start = 0
    a = 0
    b = fac2 * fac1

    while a < len(res):
        ende = start + fac1       # Python slicing is exclusive, so no -1
        # Select the slice and replicate it `fac2` times
        block = np.tile(next_prices[start:ende, :], (fac2, 1))
        res[a:b, :] = block

        start = ende
        a = b
        b += fac1 * fac2

    return res
