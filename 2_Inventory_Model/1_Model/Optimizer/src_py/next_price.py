import numpy as np
from numba import njit
num_regimes = 2

# @njit
def next_price(belief_vector: np.ndarray, 
               regime_transition_probabilities: np.ndarray
               ) -> np.ndarray:
    belief_size = belief_vector.shape[0]
    result = np.zeros((belief_size, len(regime_transition_probabilities)))

    for i in range(belief_size):
        for j in range(len(regime_transition_probabilities)):
            for r in range(num_regimes):
                result[i, j] += belief_vector[i, r] * regime_transition_probabilities[j, r]

    if belief_size == 1:  # case for t=1, for consistency.
        result = result.T

    return result