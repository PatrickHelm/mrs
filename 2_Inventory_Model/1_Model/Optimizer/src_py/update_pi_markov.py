import numpy as np
from numba import njit

CEC = 1
R1 = 2
R2 = 3


@njit(cache=True)
def _update_pi_markov_core(R1_mat, R2_mat, NoP, t, PI_past_R1, PI_past_R2,
                           k, method_code, initial_pi_r1, initial_pi_r2):
    """Core Bayesian update loop for Markov model (t >= 2).

    Parameters
    ----------
    R1_mat, R2_mat : np.ndarray (NoP, NoP)
        Transition matrices.
    NoP : int
        Number of prices.
    t : int
        Current period (1-based).
    PI_past_R1, PI_past_R2 : np.ndarray
        Previous period belief matrices.
    k : np.ndarray (2, 2)
        Regime transition probabilities.
    method_code : int
        Solution method code.
    initial_pi_r1, initial_pi_r2 : float
        Initial beliefs for CEC fallback.

    Returns
    -------
    PI_markov1, PI_markov2 : np.ndarray, shape (NoP^(t-1), NoP)
    """
    result_rows = 1
    for _ in range(t - 1):
        result_rows *= NoP

    PI_markov1 = np.zeros((result_rows, NoP))
    PI_markov2 = np.zeros((result_rows, NoP))

    i = 0  # index for PI_past (0-based)
    j = 0  # index for result (0-based)

    while j < result_rows:
        for last_price in range(NoP):
            for p_t in range(NoP):
                if method_code == CEC:
                    PI_markov1[j, p_t] = initial_pi_r1
                    PI_markov2[j, p_t] = initial_pi_r2
                elif method_code == R1:
                    PI_markov1[j, p_t] = 1.0
                    PI_markov2[j, p_t] = 0.0
                elif method_code == R2:
                    PI_markov1[j, p_t] = 0.0
                    PI_markov2[j, p_t] = 1.0
                else:
                    # Bayesian update
                    denom = (PI_past_R1[i, last_price] * R1_mat[last_price, p_t]
                             + PI_past_R2[i, last_price] * R2_mat[last_price, p_t])

                    if denom == 0.0:
                        PI_markov1[j, p_t] = PI_past_R1[i, last_price]
                        PI_markov2[j, p_t] = PI_past_R2[i, last_price]
                    else:
                        num_r1 = (PI_past_R1[i, last_price] * k[0, 0] * R1_mat[last_price, p_t]
                                  + PI_past_R2[i, last_price] * k[0, 1] * R2_mat[last_price, p_t])
                        num_r2 = (PI_past_R1[i, last_price] * k[1, 0] * R1_mat[last_price, p_t]
                                  + PI_past_R2[i, last_price] * k[1, 1] * R2_mat[last_price, p_t])
                        PI_markov1[j, p_t] = num_r1 / denom
                        PI_markov2[j, p_t] = num_r2 / denom

                        # NaN check (already handled by denom==0 above, but keep for safety)
                        if np.isnan(PI_markov1[j, p_t]):
                            PI_markov1[j, p_t] = PI_past_R1[i, last_price]
                        if np.isnan(PI_markov2[j, p_t]):
                            PI_markov2[j, p_t] = PI_past_R2[i, last_price]
            j += 1
        i += 1

    return PI_markov1, PI_markov2


def update_pi_markov(R1_mat, R2_mat, prices, t, PI_markov_list, k, method_code):
    """Compute Markov belief update PI_markov for period t.

    Parameters
    ----------
    R1_mat, R2_mat : np.ndarray (NoP, NoP)
        Markov transition matrices.
    prices : np.ndarray
        Vector of possible prices.
    t : int
        Current period (1-based).
    PI_markov_list : list or np.ndarray
        If t==1: np.ndarray of shape (2,) — initial beliefs [pi_R1, pi_R2]
        If t>=2: list of [PI_R1, PI_R2] pairs for each past period
    k : np.ndarray (2, 2)
        Regime transition probabilities.
    method_code : int
        Solution method code.

    Returns
    -------
    list : [PI_markov1, PI_markov2]
        For t==1: [scalar, scalar]
        For t>=2: [np.ndarray, np.ndarray] of shape (NoP^(t-1), NoP)
    """
    NoP = len(prices)

    if t == 1:
        # For t=1, belief values are given directly
        return [PI_markov_list[0], PI_markov_list[1]]

    # Get past beliefs
    PI_past_R1 = PI_markov_list[t - 2][0]  # beliefs in t-1 for R1
    PI_past_R2 = PI_markov_list[t - 2][1]  # beliefs in t-1 for R2

    # For t==2, PI_past is scalar -> replicate to (1, NoP) matrix
    if t == 2:
        PI_past_R1 = np.full((1, NoP), PI_past_R1)
        PI_past_R2 = np.full((1, NoP), PI_past_R2)

    # Get initial beliefs for CEC fallback
    initial_pi_r1 = float(PI_markov_list[0][0])
    initial_pi_r2 = float(PI_markov_list[0][1])

    PI_m1, PI_m2 = _update_pi_markov_core(
        R1_mat, R2_mat, NoP, t, PI_past_R1, PI_past_R2,
        k, method_code, initial_pi_r1, initial_pi_r2
    )

    return [PI_m1, PI_m2]
