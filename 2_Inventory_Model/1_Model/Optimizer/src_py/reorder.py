import numpy as np
from numba import njit


@njit(cache=True)
def aux_matrix(next_prices, prices_len, I_max, t):
    """Create expanded matrix with copies of next_prices to fit observed states.

    Equivalent to MATLAB aux_matrix.m (non-Markov version).

    Parameters
    ----------
    next_prices : np.ndarray, shape (n_filtered, n_prices)
        Filtered next-period price probabilities.
    prices_len : int
        Number of different prices (numel(prices)).
    I_max : int
        Maximum inventory level.
    t : int
        Current time period (1-based).

    Returns
    -------
    res : np.ndarray, shape (prices_len^(t-1) * (I_max+1), prices_len)
    """
    I_max_adj = I_max + 1  # account for 0-based inventory levels
    fac1 = prices_len
    fac2 = I_max_adj * fac1

    # Calculate result size
    result_rows = 1
    for _ in range(t - 1):
        result_rows *= fac1
    result_rows *= I_max_adj

    res = np.zeros((result_rows, fac1))
    start = 0
    a = 0
    b = fac2

    while a < result_rows:
        end_idx = start + fac1  # rows to copy: start:end_idx
        # Replicate the fac1 rows I_max_adj times
        for rep in range(I_max_adj):
            offset = a + rep * fac1
            for row in range(fac1):
                for col in range(fac1):
                    res[offset + row, col] = next_prices[start + row, col]
        start = end_idx
        a = b
        b = b + fac2

    return res
