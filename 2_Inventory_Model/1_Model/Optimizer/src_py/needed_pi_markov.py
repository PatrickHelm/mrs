import numpy as np
from itertools import product

def needed_pi_markov(PI_markov, observed_price, last_price, prices):
    t_max = len(PI_markov)
    PI_needed = [None] * t_max

    past_price = (prices == last_price)
    actual_price = (prices == observed_price)

    for t in range(t_max):
        R1 = PI_markov[t][0]
        R2 = PI_markov[t][1]
        if t == 0:
            PI_needed[t] = np.column_stack((R1[past_price, actual_price], R2[past_price, actual_price]))
        else:
            counter = prices
            for _ in range(t):
                counter = np.array(list(product(counter, prices)))
                counter = counter[:, ::-1]  # flip columns
            ind = (counter[:, 0] == last_price) & (counter[:, 1] == observed_price)
            PI_neededR1 = R1[ind, :]
            PI_neededR2 = R2[ind, :]
            PI_neededR1 = PI_neededR1.reshape(-1, 1)
            PI_neededR2 = PI_neededR2.reshape(-1, 1)
            PI_needed[t] = np.hstack((PI_neededR1, PI_neededR2))

    return PI_needed