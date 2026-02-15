import numpy as np
from numba import njit


@njit(cache=True)
def _next_price_markov_t1(R1, R2, NoP, pi_r1, pi_r2):
    """Compute next price probabilities for t=1."""
    nP = np.zeros((NoP, NoP))
    for p_t in range(NoP):
        for p_future in range(NoP):
            nP[p_t, p_future] = pi_r1 * R1[p_t, p_future] + pi_r2 * R2[p_t, p_future]
    return nP


@njit(cache=True)
def _next_price_markov_core(R1, R2, NoP, t, PI_R1, PI_R2):
    """Compute next price probabilities for t >= 2.

    Parameters
    ----------
    R1, R2 : np.ndarray (NoP, NoP)
        Transition matrices.
    NoP : int
        Number of prices.
    t : int
        Current period (1-based in MATLAB; here we use the same value).
    PI_R1, PI_R2 : np.ndarray
        Belief matrices for R1 and R2 at time t.
    """
    total_rows = 1
    for _ in range(t):
        total_rows *= NoP
    nP = np.zeros((total_rows, NoP))

    j = 0  # index for resulting matrix (0-based)
    i = 0  # index for PI_R1 and PI_R2 (0-based)
    while j < total_rows:
        for p_t in range(NoP):
            for p_future in range(NoP):
                nP[j, p_future] = (PI_R1[i, p_t] * R1[p_t, p_future]
                                   + PI_R2[i, p_t] * R2[p_t, p_future])
            j += 1
        i += 1
        if i % NoP == 0:
            i = 0

    return nP


def next_price_markov(R1, R2, prices, t, PI_markov):
    """Calculate probabilities for all possible prices in t+1 (Markov version).

    Parameters
    ----------
    R1, R2 : np.ndarray (NoP, NoP)
        Markov transition matrices for regime 1 and 2.
    prices : np.ndarray
        Vector of possible prices.
    t : int
        Current time period (1-based, matching MATLAB convention).
    PI_markov : varies
        If t==1: np.ndarray of shape (2,) — [pi_R1, pi_R2]
        If t>=2: list of lists, PI_markov[t-1] = [PI_R1_matrix, PI_R2_matrix]

    Returns
    -------
    nP_markov : np.ndarray, shape (NoP^t, NoP)
    """
    NoP = len(prices)

    if t == 1:
        return _next_price_markov_t1(R1, R2, NoP, PI_markov[0], PI_markov[1])
    else:
        PI_R1 = PI_markov[t - 1][0]  # belief at time t for R1
        PI_R2 = PI_markov[t - 1][1]  # belief at time t for R2
        return _next_price_markov_core(R1, R2, NoP, t, PI_R1, PI_R2)
