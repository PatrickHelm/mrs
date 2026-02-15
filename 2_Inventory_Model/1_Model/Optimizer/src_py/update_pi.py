import numpy as np
from numba import njit

# Solution method codes
EXACT = 0
CEC = 1
R1 = 2
R2 = 3
OLFC = 4


@njit(cache=True)
def update_pi(k, m, pi, P_pr_reg, t, method_code):
    """Bayesian belief update: compute PI_{t} from PI_{t-1}.

    Parameters
    ----------
    k : np.ndarray, shape (2, m)
        Regime transition probability matrix.
    m : int
        Number of regimes.
    pi : np.ndarray, shape (n_prev, m)
        Beliefs from previous period (PI{t-1}).
    P_pr_reg : np.ndarray, shape (n_prices, 2)
        Price probability distributions per regime.
    t : int
        Current time period (1-based).
    method_code : int
        0=exact, 1=CEC, 2=R1, 3=R2.

    Returns
    -------
    res : np.ndarray, shape (n_prices^(t-1), m)
    """
    n_prev = pi.shape[0]
    n_prices = P_pr_reg.shape[0]

    # Size of result
    result_rows = 1
    for _ in range(t - 1):
        result_rows *= n_prices

    int_res = np.zeros((result_rows, m))
    index = 0

    for l in range(n_prev):
        for i in range(n_prices):
            for j in range(m):
                if method_code == CEC:
                    int_res[index, j] = pi[l, j]
                elif method_code == R1:
                    if j == 0:
                        int_res[index, 0] = 1.0
                    # j==1 stays 0
                elif method_code == R2:
                    if j == 1:
                        int_res[index, 1] = 1.0
                    # j==0 stays 0
                else:
                    # Bayesian update (exact / OLFC)
                    pi_zaehler = 0.0
                    pi_nenner = 0.0
                    for n in range(2):
                        pi_zaehler += pi[l, n] * k[n, j] * P_pr_reg[i, n]
                        pi_nenner += pi[l, n] * P_pr_reg[i, n]
                    int_res[index, j] = pi_zaehler / pi_nenner
            index += 1

    return int_res
