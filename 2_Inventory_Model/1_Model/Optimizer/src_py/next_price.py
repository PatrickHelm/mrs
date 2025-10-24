import numpy as np

def next_price(PI, P_pr_reg):
    Pi_size = PI.shape
    res = np.zeros((Pi_size[0], len(P_pr_reg)))

    for i in range(Pi_size[0]):
        for j in range(len(P_pr_reg)):
            res[i, j] = PI[i, 0] * P_pr_reg[j][0] + PI[i, 1] * P_pr_reg[j][1]

    if Pi_size[0] == 1:  # case for t=1, for consistency.
        res = res.T

    return res