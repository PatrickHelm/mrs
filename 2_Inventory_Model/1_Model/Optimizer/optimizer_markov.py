"""
Optimizer_Markov.py
Translated from Optimizer_Markov.m (MATLAB) -> Python

Notes:
- Requires: numpy, pandas
- The heavy function `backward_recursion_markov` is provided as a placeholder/stub and must be replaced
  with a faithful Python implementation of the original MATLAB function.
- If the Excel file 'Stochastic_Regimes_MR-MO.xlsx' is missing, the script falls back to default placeholders.
"""

import time
import numpy as np
import pandas as pd
import os
from src_py.backward_recursion_markov import backward_recursion_markov


def main():
    # _ _SOLUTION METHOD _ _
    solution_method = 'R1'  # 'exact', 'CEC', 'R1', 'R2'

    # _ _SUBOPTIMAL DECISIONS _ _
    excel_path = '2_Inventory_Model/1_Model/Optimizer/Stochastic_Regimes_MR-MO.xlsx'
    df_excel = pd.read_excel(excel_path, sheet_name=0, header=None)

    # The MATLAB ranges:
    # suboptimal_decision_vector_R1 = xlsread(...,'B4:B8')  -> B4..B8 is rows 4..8 of column B
    # suboptimal_decision_vector_R2 = xlsread(...,'C4:C8')
    # suboptimal_decision_matrix_CEC = xlsread(...,'F4:J8') -> rows 4..8, columns F..J (5 columns)
    if df_excel is not None:
        # pandas read with header=None -> first row is index 0 corresponding to Excel row 1
        # So Excel row r maps to df_excel row r-1
            # Extract ranges using 0-based indices
            r_start = 4 - 1
            r_end = 8 - 1
            # Column B is index 1, C is index 2, F..J is 5..9
            suboptimal_decision_vector_R1 = df_excel.iloc[r_start:r_end+1, 1].astype(float).to_numpy()
            suboptimal_decision_vector_R2 = df_excel.iloc[r_start:r_end+1, 2].astype(float).to_numpy()
            suboptimal_decision_matrix_CEC = df_excel.iloc[r_start:r_end+1, 5:10].astype(float).to_numpy()

    # _ _INVENTORY COST _ _
    variables = {}
    variables['t_max'] = 4
    variables['ch'] = 1
    variables['cp'] = 40

    # _ _DEMAND _ _
    variables['d_max'] = 30
    variables['I_max'] = 120
    demand_dist = 'uniform'  # 'uniform', 'normal', 'poisson', 'nbin', 'manual'
    variables['normal'] = {'mu': 15, 'sigma': 3}
    variables['poisson'] = {'lambda': 1}
    variables['nbin'] = {'r': 22.5, 'p': 0.6}
    start_inventory = 0

    # _ _PRICE PROCESS _ _
    variables['k'] = np.array([[0.5, 0.5], [0.5, 0.5]])
    variables['m'] = 2
    beliefs = [0.1, 0.3, 0.5, 0.7, 0.9]
    variables['prices'] = [10, 15, 20, 25, 30]
    price_dist = ['MR', 'MO']
    variables['P_pr_reg'] = 0  # defined for compatibility

    # Emulate MATLAB diary by writing prints to console (user may redirect)
    # Outer loops: MATLAB loops over i=1:numel(variables.prices), h and j loops are commented to run only first element (1:1)
    num_prices = len(variables['prices'])
    for i in range(num_prices):
        observed_price = variables['prices'][i]
        for h in range(1):  # MATLAB uses 1:1 in the provided script; change to range(num_prices) if you want full sweep
            last_price = variables['prices'][h]
            for j in range(1):  # MATLAB uses 1:1; change to range(len(beliefs)) to iterate all beliefs
                pi1 = beliefs[j]
                variables['pi'] = [pi1, 1.0 - pi1]

                allSave = [None] * variables['t_max']

                # MATLAB loop: for T = variables.t_max
                for T in [variables['t_max']]:
                    variables['t_max'] = T
                    t0 = time.time()
                    periods = backward_recursion_markov(
                        variables, observed_price, last_price, start_inventory,
                        demand_dist, price_dist, solution_method,
                        suboptimal_decision_matrix_CEC, suboptimal_decision_vector_R1, suboptimal_decision_vector_R2
                    )
                    if periods:
                        # periods expected as list where index 0 behaves like periods{1} in MATLAB
                        p0 = periods[0]
                        print(f"<------------Planning horizon : {T} periods --------------------->")
                        print("-" * 70)
                        print(f"For planning horizon of {T} periods, and observed price sequence: {last_price} --> {observed_price} and initial belief {variables['pi']}:")
                        print(f"Optimal Order Decision: {p0.get('period_min_order')}")
                        print(f"Minimum Total Expected Cost: {p0.get('period_min_cost')}")
                        print("-" * 70)
                        # Safe gap calculations (avoid division by zero)
                        min_cost = p0.get('period_min_cost', 0.0) or 0.0

                        cec_cost = p0.get('suboptimal_cost_CEC', float('nan'))
                        if min_cost != 0:
                            gap_cec = (cec_cost / min_cost) * 100.0 - 100.0
                        else:
                            gap_cec = float('nan')
                        print(f"Suboptimal Cost (order decision (CEC): {suboptimal_decision_matrix_CEC}): {cec_cost}")
                        print(f"Optimality gap in % total cost (CEC): {gap_cec}")
                        print("-" * 70)

                        r1_cost = p0.get('suboptimal_cost_R1', float('nan'))
                        if min_cost != 0:
                            gap_r1 = (r1_cost / min_cost) * 100.0 - 100.0
                        else:
                            gap_r1 = float('nan')
                        print(f"Suboptimal Cost (order decision (R1): {suboptimal_decision_vector_R1}): {r1_cost}")
                        print(f"Optimality gap in % total cost (R1): {gap_r1}")
                        print("-" * 70)

                        r2_cost = p0.get('suboptimal_cost_R2', float('nan'))
                        if min_cost != 0:
                            gap_r2 = (r2_cost / min_cost) * 100.0 - 100.0
                        else:
                            gap_r2 = float('nan')
                        print(f"Suboptimal Cost (order decision (R2): {suboptimal_decision_vector_R2}): {r2_cost}")
                        print(f"Optimality gap in % total cost (R2): {gap_r2}")

                    t_elapsed = time.time() - t0
                    print(f"(Elapsed time: {t_elapsed:.4f} seconds)")
                    # MATLAB stored allSave{T} = {periods}; here we store at index T-1
                    if 1 <= T <= len(allSave):
                        allSave[T-1] = periods
    # End loops
    print("Finished run.")


if __name__ == '__main__':
    main()