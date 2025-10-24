from itertools import product
import numpy as np
from .get_demand_dist import get_demand_dist
from .get_markov_model import getMarkovModel
from .get_mins import get_mins
from .allocate import allocate
from .backward_cost_function import backward_cost_function
from .single_period_cost import single_period_cost
from .update_pi_markov import update_pi_markov
from .next_price_markov import next_price_markov
from .reorder_markov import reorder_markov
from .needed_pi_markov import needed_pi_markov
from .get_next_prices_markov import get_next_prices_markov

def backward_recursion_markov(variables, observed_price, last_price, start_inventory, dist, priceDist, solution_method, suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2):
    variables_tuple = allocate(variables)
    ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma, poisson_lambda, nbin_r, nbin_p = variables_tuple
    demand_dist = get_demand_dist(dist, d_max, variables_tuple)

    PI_markov = np.empty(t_max)
    nP_markov = np.empty(t_max)
    markovmodel = getMarkovModel(priceDist)

    for t in range(t_max + 1):
        if t == 0:
            PI_markov[t] = update_pi_markov(markovmodel, prices, t + 1, pi, k, solution_method)  # Estimated values for regime 1 or regime 2 in period 1 (given); here [0.5 0.5]
            nP_markov[t] = next_price_markov(markovmodel, prices, t + 1, pi)
        else:
            PI_markov[t] = update_pi_markov(markovmodel, prices, t + 1, PI_markov, k, solution_method)
            nP_markov[t] = next_price_markov(markovmodel, prices, t + 1, PI_markov)

    PI_markov = PI_markov[1:]
    nP_markov = nP_markov[:-1]
    PI_needed_markov = needed_pi_markov(PI_markov, observed_price, last_price, prices)

    if np.array_equal(PI_needed_markov[0], [0, 0]):
        print('No reachable state')
        return {}

    PI_needed = [pi[:, 0] for pi in PI_needed_markov]

    inv = np.arange(I_max + 1)

    periods = [None] * (t_max + 1)

    for t in range(t_max, 0, -1):  # starting in the last period
        if t == t_max:
            all_comb = np.array(list(product(prices, inv, PI_needed[t - 1])))
            period_matrix = np.zeros((len(all_comb), I_max + 1))
            period_min_cost = np.zeros(len(all_comb))
            period_min_order = np.zeros(len(all_comb), dtype=int)

            for s in range(len(all_comb)):
                for y in range(I_max + 1):
                    period_matrix[s, y] = single_period_cost(y, cp, ch, all_comb[s, 0], d_max, all_comb[s, 1], demand_dist, dist)
                period_min_cost[s] = np.min(period_matrix[s, :])
                period_min_order[s] = np.argmin(period_matrix[s, :])

            periods[t] = {
                'period_matrix': period_matrix,
                'period_min_cost': period_min_cost,
                'period_min_order': period_min_order
            }
        else:
            if t != 1:
                all_comb = np.array(list(product(prices, inv, PI_needed[t - 1])))
            else:
                all_comb = np.array(list(product([observed_price], [start_inventory], PI_needed[t - 1])))

            period_matrix = np.zeros((len(all_comb), I_max + 1))
            period_min_cost = np.zeros(len(all_comb))
            period_min_order = np.zeros(len(all_comb), dtype=int)

            next_prices = get_next_prices_markov(nP_markov[t], all_comb[:, 0], observed_price, last_price, prices, t)
            if t != 1:
                next_prices = reorder_markov(next_prices, prices, I_max, t)

            for s in range(len(all_comb)):
                mins = get_mins(PI_needed, prices, t, all_comb[s, 2], all_comb[s, 0], periods[t + 1]['period_min_cost'], inv)
                for y in range(I_max + 1):
                    period_matrix[s, y] = backward_cost_function(y, cp, ch, all_comb[s, 0], d_max, I_max, all_comb[s, 1], next_prices[s, :], mins, demand_dist, dist)
                period_min_cost[s] = np.min(period_matrix[s, :])
                period_min_order[s] = np.argmin(period_matrix[s, :])

           
    periods[t] = {
            'PI': PI_markov[t],
            'PI_needed': PI_needed[t],
            'period_matrix': period_matrix,
            'period_min_cost': period_min_cost,
            'period_min_order': period_min_order,
            'suboptimal_cost_CEC': period_matrix[-1, suboptimal_decision_CEC],
            'suboptimal_cost_R1': period_matrix[-1, suboptimal_decision_R1],
            'suboptimal_cost_R2': period_matrix[-1, suboptimal_decision_R2],
            'period_min_cost': period_min_cost,
            'period_min_order': period_min_order,
            'next_price':  nP_markov[t]
        }
 

    return periods