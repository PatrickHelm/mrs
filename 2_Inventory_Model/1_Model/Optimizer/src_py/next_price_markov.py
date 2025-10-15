import numpy as np

def next_price_markov(markovmodel, prices, t, PI_markov):
    NoP = len(prices)
    nP_markov = np.zeros((NoP**t, NoP))
    R1 = markovmodel[0]  # probability transition matrix for R1
    R2 = markovmodel[1]  # probability transition matrix for R2

    if t == 1:
        for p_t in range(NoP):  # p_t is price observation in period t
            for p_future in range(NoP):  # p_future = p_t+1
                nP_markov[p_t, p_future] = PI_markov[0] * R1[p_t, p_future] + PI_markov[1] * R2[p_t, p_future]
    else:
        PI_R1 = PI_markov[t-1][0]  # next period regime belief at time t for R1
        PI_R2 = PI_markov[t-1][1]  # next period regime belief at time t for R2
        j = 0  # index for resulting matrix
        i = 0  # index PI_R1 and PI_R2
        while j < NoP**t:
            for p_t in range(NoP):
                for p_future in range(NoP):
                    nP_markov[j, p_future] = PI_R1[i, p_t] * R1[p_t, p_future] + PI_R2[i, p_t] * R2[p_t, p_future]
                j += 1
            if (i + 1) % NoP == 0:
                i = -1
            i += 1

    return nP_markov