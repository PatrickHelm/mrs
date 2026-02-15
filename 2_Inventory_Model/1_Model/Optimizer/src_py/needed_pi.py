import numpy as np
from itertools import product as iterproduct


def _combvec_flipped(prices, num_dims):
    """Equivalent to MATLAB's iterative combvec + fliplr(var') pattern.

    itertools.product already produces standard lexicographic order,
    which matches MATLAB's combvec+fliplr result. No fliplr needed.
    Returns shape (NoP^num_dims, num_dims).
    """
    tuples = np.array(list(iterproduct(prices, repeat=num_dims)), dtype=np.float64)
    return tuples


def needed_pi(PI_list, observed_price, prices):
    """Return reachable beliefs conditioned on the observed price in period 1.

    Parameters
    ----------
    PI_list : list of np.ndarray
        PI_list[t] is the belief array for period t+1 (0-indexed).
        PI_list[0] shape (1, m) or (1,), PI_list[1] shape (n_prices, m), etc.
    observed_price : float
        Price observed in period 1.
    prices : np.ndarray
        Vector of all possible prices.

    Returns
    -------
    res : list of np.ndarray
        res[t] contains the reachable belief values (column 0 of PI) for period t+1.
    """
    t_max = len(PI_list)
    res = [None] * t_max

    # Boolean mask for observed price
    actual_price_mask = (prices == observed_price)

    for i in range(t_max):
        if i == 0:
            # Period 1: initial belief (first element of first column)
            pi_arr = PI_list[0]
            if pi_arr.ndim == 1:
                res[0] = np.array([pi_arr[0]])
            else:
                res[0] = np.array([pi_arr[0, 0]])

        elif i == 1:
            # Period 2: select row matching observed price
            pi_arr = PI_list[1]
            if pi_arr.ndim == 1:
                res[1] = np.array([pi_arr[actual_price_mask][0]])
            else:
                res[1] = pi_arr[actual_price_mask, 0]

        else:
            # Period >= 3: use combvec + fliplr pattern
            # MATLAB i_matlab = i+1, loop runs j=1:i_matlab-2, building i_matlab-1 = i dims.
            # Result: NoP^i rows × i cols, flipped. Filter on column 0 == observed_price.
            var = _combvec_flipped(prices, i)
            ind = (var[:, 0] == observed_price)
            pi_arr = PI_list[i]
            if pi_arr.ndim == 1:
                res[i] = pi_arr[ind]
            else:
                res[i] = pi_arr[ind, 0]

    return res
