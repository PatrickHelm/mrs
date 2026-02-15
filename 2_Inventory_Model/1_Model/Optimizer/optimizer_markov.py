"""Inventory optimizer — Markov version.

Translated from Optimizer_Markov.m.
"""
import time
import numpy as np
import openpyxl

from src_py.backward_recursion_markov import backward_recursion_markov

# --- Solution method codes ---
EXACT = 0
CEC = 1
R1 = 2
R2 = 3

METHOD_MAP = {'exact': EXACT, 'CEC': CEC, 'R1': R1, 'R2': R2}

# --- Demand distribution codes ---
DIST_MAP = {'uniform': 0, 'normal': 1, 'poisson': 2, 'nbin': 3, 'manual': 4}


def read_excel_range(filepath, sheet_idx, range_str):
    """Read a rectangular range from an Excel file."""
    wb = openpyxl.load_workbook(filepath, data_only=True)
    ws = wb.worksheets[sheet_idx]
    data = []
    for row in ws[range_str]:
        row_data = []
        for cell in row:
            val = cell.value
            row_data.append(float(val) if val is not None else 0.0)
        data.append(row_data)
    wb.close()
    arr = np.array(data, dtype=np.float64)
    if arr.shape[1] == 1:
        return arr.ravel()
    return arr


def main():
    # __ SOLUTION METHOD __
    solution_method = 'R1'
    method_code = METHOD_MAP[solution_method]

    # __ SUBOPTIMAL DECISIONS __
    suboptimal_decision_vector_R1 = read_excel_range(
        'Stochastic_Regimes_MR-MO.xlsx', 0, 'B4:B8')
    suboptimal_decision_vector_R2 = read_excel_range(
        'Stochastic_Regimes_MR-MO.xlsx', 0, 'C4:C8')
    suboptimal_decision_matrix_CEC = read_excel_range(
        'Stochastic_Regimes_MR-MO.xlsx', 0, 'F4:J8')

    # __ INVENTORY COST __
    t_max = 4
    ch = 1.0
    cp = 40.0

    # __ DEMAND __
    d_max = 30
    I_max = 120
    demand_dist = 'uniform'
    dist_code = DIST_MAP[demand_dist]
    normal_mu = 15.0
    normal_sigma = 3.0
    poisson_lambda = 1.0
    nbin_r = 22.5
    nbin_p = 0.6
    start_inventory = 0

    # __ PRICE PROCESS __
    k = np.array([[0.5, 0.5], [0.5, 0.5]])
    m = 2
    beliefs = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    prices = np.array([10.0, 15.0, 20.0, 25.0, 30.0])
    price_dist = ['MR', 'MO']
    P_pr_reg = 0.0  # not used for Markov case, but must be defined

    # __ BUILD VARIABLES DICT __
    variables = {
        'ch': ch, 'cp': cp, 'd_max': d_max, 'I_max': I_max,
        'k': k, 'm': m, 'pi': None,  # set per iteration
        'P_pr_reg': P_pr_reg, 't_max': t_max, 'prices': prices,
        'normal_mu': normal_mu, 'normal_sigma': normal_sigma,
        'poisson_lambda': poisson_lambda, 'nbin_r': nbin_r, 'nbin_p': nbin_p,
    }

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
    main()
