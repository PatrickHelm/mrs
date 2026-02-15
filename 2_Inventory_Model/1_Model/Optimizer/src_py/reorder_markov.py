import numpy as np
from numba import njit


@njit(cache=True)
def reorder_markov(next_prices, prices_len, I_max, t):
    """Create expanded matrix with copies of next_prices for Markov model.

    Equivalent to MATLAB reorderMarkov.m.

    Parameters
    ----------
    next_prices : np.ndarray, shape (n_filtered, n_prices)
        Filtered next-period price probabilities.
    prices_len : int
        Number of different prices.
    I_max : int
        Maximum inventory level.
    t : int
        Current time period (1-based).

    Returns
    -------
    res : np.ndarray, shape (prices_len^t * (I_max+1), prices_len)
    """
    I_max_adj = I_max + 1
    fac1 = prices_len
    fac2 = I_max_adj * fac1

    # Calculate result size
    result_rows = 1
    for _ in range(t):
        result_rows *= fac1
    result_rows *= I_max_adj

    res = np.zeros((result_rows, fac1))
    start = 0
    a = 0
    b = fac2 * fac1

    while a < result_rows:
        end_idx = start + fac1
        # Replicate the fac1 rows fac2 times (not I_max_adj like aux_matrix)
        for rep in range(fac2):
            offset = a + rep * fac1
            for row in range(fac1):
                for col in range(fac1):
                    res[offset + row, col] = next_prices[start + row, col]
        start = end_idx
        a = b
        b = b + fac1 * fac2

    return res
