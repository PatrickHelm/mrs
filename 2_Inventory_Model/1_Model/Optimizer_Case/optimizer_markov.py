"""Inventory optimizer — Case Markov version.

Translated from Optimizer_Case/Optimizer_Markov.m.
Loads all parameters from a case data directory (default: 3_Case/Zinc/).
"""
import sys
import os
import time
import numpy as np

# Add parent Optimizer directory to path so we can import src_py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Optimizer'))
from src_py.backward_recursion_markov import backward_recursion_markov
from load_case import load_case_markov


def main(case_dir=None):
    if case_dir is None:
        case_dir = os.path.join(os.path.dirname(__file__), '..', '..',
                                '3_Case', 'Zinc')

    (variables, beliefs, suboptimal_decision_vector_R1,
     suboptimal_decision_vector_R2, suboptimal_decision_matrix_CEC,
     _, method_code, dist_code, start_inventory
     ) = load_case_markov(case_dir)

    prices = variables['prices']
    t_max = variables['t_max']
    # price_dist labels not needed — R1/R2 are loaded directly from files
    price_dist = ['R1', 'R2']

    # __ CALCULATIONS START __
    for i in range(len(prices)):
        observed_price = prices[i]
        for h in range(1):  # range(len(prices))
            last_price = prices[h]
            for j in range(1):  # range(len(beliefs))
                pi1 = beliefs[j]
                variables['pi'] = np.array([pi1, 1.0 - pi1])

                suboptimal_decision_CEC = int(suboptimal_decision_matrix_CEC[j, i])
                suboptimal_decision_R1 = int(suboptimal_decision_vector_R1[i])
                suboptimal_decision_R2 = int(suboptimal_decision_vector_R2[i])

                for T in [t_max]:
                    variables['t_max'] = T
                    tic = time.time()

                    periods = backward_recursion_markov(
                        variables, observed_price, last_price, start_inventory,
                        dist_code, price_dist, method_code,
                        suboptimal_decision_CEC, suboptimal_decision_R1,
                        suboptimal_decision_R2
                    )

                    toc = time.time()

                    if periods is not None:
                        print(f'<------------Planning horizon : {T} periods --------------------->')
                        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
                        print(f'For planning horizon of {T} periods, and observed price '
                              f'sequence: {last_price}-->{observed_price} and initial belief '
                              f'{variables["pi"]}:')
                        print(f'Optimal Order Decision: {periods[0]["period_min_order"]}')
                        print(f'Minimum Total Expected Cost: {periods[0]["period_min_cost"]}')
                        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
                        print(f'Suboptimal Cost (order decision (CEC): '
                              f'{suboptimal_decision_CEC}): '
                              f'{periods[0]["suboptimal_cost_CEC"]}')
                        opt_gap_CEC = (periods[0]['suboptimal_cost_CEC']
                                       / periods[0]['period_min_cost'][0] * 100 - 100)
                        print(f'Optimality gap in % total cost (CEC): {opt_gap_CEC}')
                        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
                        print(f'Suboptimal Cost (order decision (R1): '
                              f'{suboptimal_decision_R1}): '
                              f'{periods[0]["suboptimal_cost_R1"]}')
                        opt_gap_R1 = (periods[0]['suboptimal_cost_R1']
                                      / periods[0]['period_min_cost'][0] * 100 - 100)
                        print(f'Optimality gap in % total cost (R1): {opt_gap_R1}')
                        print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
                        print(f'Suboptimal Cost (order decision (R2): '
                              f'{suboptimal_decision_R2}): '
                              f'{periods[0]["suboptimal_cost_R2"]}')
                        opt_gap_R2 = (periods[0]['suboptimal_cost_R2']
                                      / periods[0]['period_min_cost'][0] * 100 - 100)
                        print(f'Optimality gap in % total cost (R2): {opt_gap_R2}')

                    print(f'Elapsed time: {toc - tic:.3f} seconds')


if __name__ == '__main__':
    case = sys.argv[1] if len(sys.argv) > 1 else None
    main(case)
