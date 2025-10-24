
import numpy as np

def update_pi_markov(markovmodel, prices, t, PI, k, solution_method):
    NoP = len(prices)  # Number of Prices
    PI_markov = [None, None]
    R1 = markovmodel[0]  # probability transition matrix for R1
    R2 = markovmodel[1]  # probability transition matrix for R2

    if t == 1:  # for t=1 the belief update values are given
        PI_markov[0] = PI[0]
        PI_markov[1] = PI[1]
        return PI_markov

    PI_markov1 = np.zeros((NoP**(t-1), NoP))
    PI_markov2 = np.zeros((NoP**(t-1), NoP))
    PI_past_R1 = PI[t-2][0]  # belief updates in t-1 for R1
    PI_past_R2 = PI[t-2][1]  # belief updates in t-1 for R2

    if t == 2:  # blow up the size of the belief in period 1 by copying its value NoP times
        PI_past_R1 = np.tile(PI_past_R1, (1, NoP))
        PI_past_R2 = np.tile(PI_past_R2, (1, NoP))

    i = 0
    j = 0
    while j < NoP**(t-1):
        for last_price in range(NoP):  # p_past is the price in the last period p_(t-1)
            for p_t in range(NoP):  # p_t is the price in the actual period t
                if solution_method == 'CEC':
                    PI_markov1[j, p_t] = PI[0][0]
                    PI_markov2[j, p_t] = PI[0][1]
                elif solution_method == 'R1':
                    PI_markov1[j, p_t] = 1
                    PI_markov2[j, p_t] = 0
                elif solution_method == 'R2':
                    PI_markov1[j, p_t] = 0
                    PI_markov2[j, p_t] = 1
                else:
                    numerator_R1 = (PI_past_R1[i, last_price] * k[0, 0] * R1[last_price, p_t] +
                                    PI_past_R2[i, last_price] * k[0, 1] * R2[last_price, p_t])
                    denominator = (PI_past_R1[i, last_price] * R1[last_price, p_t] +
                                   PI_past_R2[i, last_price] * R2[last_price, p_t])
                    numerator_R2 = (PI_past_R1[i, last_price] * k[1, 0] * R1[last_price, p_t] +
                                    PI_past_R2[i, last_price] * k[1, 1] * R2[last_price, p_t])

                    PI_markov1[j, p_t] = numerator_R1 / denominator if denominator != 0 else PI_past_R1[i, last_price]
                    PI_markov2[j, p_t] = numerator_R2 / denominator if denominator != 0 else PI_past_R2[i, last_price]

        j += 1
        i += 1

    PI_markov = [PI_markov1, PI_markov2]
    return PI_markov