import numpy as np


def get_markov_model(price_dist):
    """Create transition probability matrices for R1 and R2.

    According to Gavirneni (2004). Each regime can be:
    - 'RW' : Random Walk
    - 'MR' : Mean-Reverting
    - 'MO' : Momentum

    Parameters
    ----------
    price_dist : list of str, length 2
        e.g. ['MR', 'MO']

    Returns
    -------
    R1 : np.ndarray (5, 5)
    R2 : np.ndarray (5, 5)
    """
    matrices = {
        'RW': np.array([
            [0.5, 0.5, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.5],
            [0.0, 0.0, 0.0, 0.5, 0.5],
        ]),
        'MR': np.array([
            [0.2, 0.8, 0.0, 0.0, 0.0],
            [0.1, 0.3, 0.6, 0.0, 0.0],
            [0.0, 0.1, 0.8, 0.1, 0.0],
            [0.0, 0.0, 0.6, 0.3, 0.1],
            [0.0, 0.0, 0.0, 0.8, 0.2],
        ]),
        'MO': np.array([
            [0.8, 0.2, 0.0, 0.0, 0.0],
            [0.6, 0.3, 0.1, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.1, 0.3, 0.6],
            [0.0, 0.0, 0.0, 0.2, 0.8],
        ]),
    }

    R1 = matrices[price_dist[0]].copy()
    R2 = matrices[price_dist[1]].copy()
    return R1, R2
