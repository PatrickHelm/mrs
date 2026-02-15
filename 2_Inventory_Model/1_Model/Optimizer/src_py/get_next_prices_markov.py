import numpy as np
from itertools import product as iterproduct


def _combvec_flipped(prices, num_dims):
    """Equivalent to MATLAB's iterative combvec + fliplr(var') pattern.

    itertools.product already produces standard lexicographic order,
    which matches MATLAB's combvec+fliplr result. No fliplr needed.
    """
    tuples = np.array(list(iterproduct(prices, repeat=num_dims)), dtype=np.float64)
    return tuples


def get_next_prices_markov(nP, p, observed_price, last_price, prices, t):
    """Compute next-period price probabilities for the Markov model.

    Filters on both last_price and observed_price.

    Parameters
    ----------
    nP : np.ndarray
        Full next-price probability matrix from next_price_markov() for period t+1.
    p : np.ndarray or float
        Current price(s) — used only for t<1 (edge case).
    observed_price : float
        Price observed in current period.
    last_price : float
        Price observed in previous period.
    prices : np.ndarray
        Vector of possible prices.
    t : int
        Current time period (1-based).

    Returns
    -------
    res : np.ndarray
        Filtered price probabilities.
    """
    if t >= 1:
        # Build all (t+1)-tuples of prices
        # MATLAB: var=prices; for i=1:t, var=combvec(var,prices); end
        # After t iterations: t+1 dimensions
        var = _combvec_flipped(prices, t + 1)
        # Filter: column 0 == last_price AND column 1 == observed_price
        pInd1 = (var[:, 0] == last_price) & (var[:, 1] == observed_price)
        res = nP[pInd1, :]
    else:
        # t < 1 edge case (shouldn't normally happen)
        for i in range(len(prices)):
            if last_price == prices[i]:
                res = nP[i:i+1, :]
                break

    return res
