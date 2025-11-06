import numpy as np
from itertools import product
from get_demand_dist import get_demand_dist
from get_mins import get_mins
from allocate import allocate
from backward_cost_function import backward_cost_function
from single_period_cost import single_period_cost
from update_pi import update_pi
from needed_pi import needed_pi
from get_next_prices import get_next_prices
from reorder import reorder
from next_price import next_price

def backward_recursion(variables, observed_price, start_inventory, dist, solution_method, suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2):
    ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma, poisson_lambda, nbin_r, nbin_p = allocate(variables)
    demand_dist = get_demand_dist(dist, d_max, variables)

    PI = [None] * t_max  # variable that will contain the belief values for t=1, 2,...t_max
    nP = [None] * t_max  # variable that will contain the next prices for t=1, 2,...t_max
    for t in range(t_max):
        if t == 0:
            PI[t] = pi  # initial beliefs for regime 1 and 2 in period 1 (given)
        else:
            if solution_method == 'OLFC':
                PI[t] = update_pi(k, m, pi, P_pr_reg, t + 1, solution_method)
            else:
                PI[t] = update_pi(k, m, PI[t - 1], P_pr_reg, t + 1, solution_method)
        nP[t] = next_price(PI[t], P_pr_reg)  # calculate probabilities for all possible prices in period t+1

    PI_needed = needed_pi(PI, observed_price, prices)  # only consider beliefs that can be reached, given the observed price in period 1

    inv = np.arange(I_max + 1)

    periods = [None] * t_max

    for t in range(t_max - 1, -1, -1):  # starting in the last period
        if t == t_max - 1:
            # only states that can be reached are under examination
            all_comb = np.array(list(product(prices, inv, PI_needed[t].T)))
            period_matrix = np.zeros((len(all_comb), I_max + 1))
            period_min_cost = np.zeros(len(all_comb))
            period_min_order = np.zeros(len(all_comb), dtype=int)
            for s in range(len(all_comb)):
                for y in range(I_max + 1):
                    period_matrix[s, y] = single_period_cost(y, cp, ch, all_comb[s, 0], d_max, all_comb[s, 1], demand_dist, dist)
                period_min_cost[s] = np.min(period_matrix[s, :])
                period_min_order[s] = np.argmin(period_matrix[s, :])
        else:
            if t != 0:
                all_comb = np.array(list(product(prices, inv, PI_needed[t].T)))
            else:
                all_comb = np.array(list(product([observed_price], [start_inventory], PI_needed[t].T)))
            all_comb_size = all_comb.shape
            period_matrix = np.zeros((all_comb_size[0], I_max + 1))
            period_min_cost = np.zeros(all_comb_size[0])
            period_min_order = np.zeros(all_comb_size[0], dtype=int)
            next_prices = get_next_prices(nP[t + 1], all_comb[:, 0], observed_price, prices, t + 1)
            if t != 0:
                next_prices = reorder(next_prices, prices, I_max, t + 1)
            for s in range(all_comb_size[0]):
                mins = get_mins(PI_needed, prices, t + 1, all_comb[s, 2], all_comb[s, 0], periods[t + 1]['period_min_cost'], inv)
                for y in range(I_max + 1):
                    period_matrix[s, y] = backward_cost_function(y, cp, ch, all_comb[s, 0], d_max, I_max, all_comb[s, 1], next_prices[s, :], mins, demand_dist, dist)
                period_min_cost[s] = np.min(period_matrix[s, :])
                period_min_order[s] = np.argmin(period_matrix[s, :])

        periods[t] = {
            'PI': PI[t],
            'PI_needed': PI_needed[t],
            'period_matrix': period_matrix,
            'period_min_cost': period_min_cost,
            'period_min_order': period_min_order,
            'suboptimal_cost_CEC': period_matrix[-1, suboptimal_decision_CEC],
            'suboptimal_cost_R1': period_matrix[-1, suboptimal_decision_R1],
            'suboptimal_cost_R2': period_matrix[-1, suboptimal_decision_R2],
'period_min_cost': period_min_cost,
'period_min_order': period_min_order,
'next_price':  nP[t]
        }

    return periods
