import numpy as np
from itertools import product as iterproduct


def _combvec_flipped(prices, num_dims):
    """Equivalent to MATLAB's iterative combvec + fliplr(var') pattern.

    itertools.product already produces standard lexicographic order,
    which matches MATLAB's combvec+fliplr result. No fliplr needed.
    """
    tuples = np.array(list(iterproduct(prices, repeat=num_dims)), dtype=np.float64)
    return tuples


def needed_pi_markov(PI_markov_list, observed_price, last_price, prices):
    """Return reachable Markov beliefs given observed price transitions.

    Parameters
    ----------
    PI_markov_list : list of [R1_matrix, R2_matrix]
        PI_markov_list[t] = [PI_R1, PI_R2] for period t+1 (0-indexed).
    observed_price : float
        Currently observed price.
    last_price : float
        Previously observed price.
    prices : np.ndarray
        Vector of possible prices.

    Returns
    -------
    PI_needed : list of np.ndarray
        PI_needed[t] has shape (n, 2) with [pi_R1, pi_R2] columns.
    """
    t_max = len(PI_markov_list)
    PI_needed = [None] * t_max

    past_price_idx = np.where(prices == last_price)[0][0]
    actual_price_idx = np.where(prices == observed_price)[0][0]

    for t in range(t_max):
        R1 = PI_markov_list[t][0]
        R2 = PI_markov_list[t][1]

        if t == 0:
            # Period 1: scalar beliefs from R1/R2 at (past_price, actual_price)
            if np.isscalar(R1):
                r1_val = R1
                r2_val = R2
            else:
                r1_val = R1[past_price_idx, actual_price_idx]
                r2_val = R2[past_price_idx, actual_price_idx]
            PI_needed[t] = np.array([[r1_val, r2_val]])
        else:
            # MATLAB loop runs j=1:t-1, building t dimensions after fliplr.
            # NoP^t rows × t cols.
            var = _combvec_flipped(prices, t + 1)
            # Filter: column 0 == last_price AND column 1 == observed_price
            ind = (var[:, 0] == last_price) & (var[:, 1] == observed_price)

            PI_needed_R1 = PI_markov_list[t][0][ind, :]
            PI_needed_R2 = PI_markov_list[t][1][ind, :]

            # Reshape to column vectors and concatenate
            PI_needed_R1_flat = PI_needed_R1.reshape(-1, 1)
            PI_needed_R2_flat = PI_needed_R2.reshape(-1, 1)
            PI_needed[t] = np.hstack([PI_needed_R1_flat, PI_needed_R2_flat])

    # Check if no reachable state
    if PI_needed[0].shape[0] == 1 and np.all(PI_needed[0] == 0):
        return PI_needed  # caller checks for [0, 0]

    return PI_needed
