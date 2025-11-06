import numpy as np
from numba import njit

# @njit
def next_price_markov(markov_model: np.ndarray, num_prices: int, t: int, belief_markov: np.ndarray):
    num_regimes = belief_markov.shape[1]
    next_price_markov = np.zeros((num_prices**t, num_prices))
    i = 0
    for j in range(num_prices**t):
        for p_t in range(num_prices):
            for p_future in range(num_prices):
                for r in range(num_regimes):
                    next_price_markov[j, p_future] += belief_markov[t-1, r, i, p_t] * markov_model[r, p_t, p_future]
        if (i + 1) % num_prices == 0:
            i = -1
        i += 1

    return next_price_markov