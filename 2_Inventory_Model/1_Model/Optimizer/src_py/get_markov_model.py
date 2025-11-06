
import numpy as np
from numba import njit

# @njit
def get_markov_model(transition_name_list: list[str], num_regimes: int, num_prices: int):
    # GETMARKOVMODEL Creates the transition probability matrix for R1 and R2 
    # MarkovModel = getMarkovModel(pDist) creates a list of length 2 containing the 
    # transition probability matrix for R1 and R2 according to Gavirneni(2004).
    # pDist itself is a list of length 2 containing the name of the model for
    # exchange rate fluctuations. This can either be
    # - 'random_walk' (random walk)
    # - 'mean_reverting' (mean-reverting)
    # - 'momentum' (momentum)

    MarkovModel = np.empty((num_regimes, num_prices, num_prices))
    for R in range(num_regimes):
        if transition_name_list[R] == 'random_walk':
            MarkovModel[R] = np.array([
                [0.5, 0.5, 0, 0, 0],
                [0.5, 0, 0.5, 0, 0],
                [0, 0.5, 0, 0.5, 0],
                [0, 0, 0.5, 0, 0.5],
                [0, 0, 0, 0.5, 0.5]
            ])
        elif transition_name_list[R] == 'mean_reverting':
            MarkovModel[R] = np.asarray([
                [0.2, 0.8, 0, 0, 0],
                [0.1, 0.3, 0.6, 0, 0],
                [0, 0.1, 0.8, 0.1, 0],
                [0, 0, 0.6, 0.3, 0.1],
                [0, 0, 0, 0.8, 0.2]
            ])
        elif transition_name_list[R] == 'momentum':
            MarkovModel[R] = np.asarray([
                [0.8, 0.2, 0, 0, 0],
                [0.6, 0.3, 0.1, 0, 0],
                [0, 0.5, 0, 0.5, 0],
                [0, 0, 0.1, 0.3, 0.6],
                [0, 0, 0, 0.2, 0.8]
            ])
    return MarkovModel
