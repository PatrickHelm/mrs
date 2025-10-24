import numpy as np

def update_pi(k, m, pi, P_pr_reg, t, solution_method):
    Pi_size = pi.shape  # size of PI{t-1}
    int_res = np.zeros((len(P_pr_reg)**(t-1), m))  # size of the result, hence of PI{t}
    index = 0
    for l in range(Pi_size[0]):
        for i in range(len(P_pr_reg)):
            for j in range(m):  # regimes
                if solution_method == 'CEC':
                    int_res[index, j] = pi[l, j]
                elif solution_method == 'R1':
                    int_res[index, 0] = 1
                elif solution_method == 'R2':
                    int_res[index, 1] = 1
                else:
                    pi_zaehler = 0
                    pi_nenner = 0
                    for n in range(2):
                        pi_zaehler += pi[l, n] * k[n, j] * P_pr_reg[i][n]
                        pi_nenner += pi[l, n] * P_pr_reg[i][n]
                    int_res[index, j] = pi_zaehler / pi_nenner
            index += 1

    res = int_res
    return res
