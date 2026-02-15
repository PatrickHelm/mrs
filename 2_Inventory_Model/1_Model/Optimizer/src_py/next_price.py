import numpy as np
from numba import njit


@njit(cache=True)
def next_price(PI, P_pr_reg):
    """Calculate probabilities for all possible prices in t+1.

    Parameters
    ----------
    PI : np.ndarray, shape (n_beliefs, 2)
        Belief matrix. Each row is [prob_regime1, prob_regime2].
    P_pr_reg : np.ndarray, shape (n_prices, 2)
        Price probability distributions per regime.

    Returns
    -------
    res : np.ndarray
        If n_beliefs == 1: shape (n_prices, 1) — column vector for consistency.
        Otherwise: shape (n_beliefs, n_prices).
    """
    n_beliefs = PI.shape[0]
    n_prices = P_pr_reg.shape[0]
    res = np.zeros((n_beliefs, n_prices))

    for i in range(n_beliefs):
        for j in range(n_prices):
            res[i, j] = PI[i, 0] * P_pr_reg[j, 0] + PI[i, 1] * P_pr_reg[j, 1]

    if n_beliefs == 1:
        # Transpose for consistency (MATLAB returns column vector for t=1)
        return res.T
    return res
