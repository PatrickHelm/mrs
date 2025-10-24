import numpy as np
from itertools import product

def get_mins(PI_needed, prices, t, act_pi, act_price, all_mins, inv):
    if t == 1:
        mins = all_mins
    else:
        block = len(prices)  # number of possible prices
        index = np.where(np.array(PI_needed[t]) == act_pi)[0][0]  # index of act_pi in PI_needed[t]
        ende = (index + 1) * block
        start = ende - block
        future_PI = PI_needed[t + 1][start:ende]  # belief updates in t+1 from act_pi
        index_2 = np.where(prices == act_price)[0][0]  # index of act_price in prices
        PI = future_PI[index_2]  # belief update in t+1 observing act_price

        # Create all combinations of prices, inv, and PI_needed[t+1]
        all_comb_future = np.array(list(product(prices, inv, PI_needed[t + 1])))
        index_3 = all_comb_future[:, 2] == PI  # logical array where belief update equals PI
        mins = all_mins[index_3]

    return mins
