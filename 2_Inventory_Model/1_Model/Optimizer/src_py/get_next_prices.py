import numpy as np
from itertools import product as iterproduct


def _combvec_flipped(prices, num_dims):
    """Equivalent to MATLAB's iterative combvec + fliplr(var') pattern.

    itertools.product already produces standard lexicographic order,
    which matches MATLAB's combvec+fliplr result. No fliplr needed.
    """
    tuples = np.array(list(iterproduct(prices, repeat=num_dims)), dtype=np.float64)
    return tuples


def get_next_prices(nP, p, observed_price, prices, t):
    """Compute next-period price probabilities conditioned on observed price.

    Parameters
    ----------
    nP : np.ndarray
        Full next-price probability matrix from next_price() for period t+1.
    p : np.ndarray or float
        Current price(s) — used only for t==1.
    observed_price : float
        Price observed in period 1.
    prices : np.ndarray
        Vector of possible prices.
    t : int
        Current time period (1-based).

    Returns
    -------
    res : np.ndarray
        Filtered price probabilities. For t==1: single row.
        For t>=2: subset of rows matching observed_price filter.
    """
    if t >= 2:
        # Build all t-tuples of prices (combvec iterated t-1 times from prices)
        # MATLAB: var=prices; for i=1:t-1, var=combvec(var,prices); end
        # After t-1 iterations: t dimensions, NoP^t rows after transpose
        var = _combvec_flipped(prices, t)
        # Filter: column 0 == observed_price
        pInd1 = (var[:, 0] == observed_price)
        res = nP[pInd1, :]
    else:
        # t == 1: find the row in nP matching p
        if np.isscalar(p):
            p_val = p
        else:
            p_val = p[0] if len(p) == 1 else p[0]

        for i in range(len(prices)):
            if p_val == prices[i]:
                res = nP[i:i+1, :]
                break

    return res
