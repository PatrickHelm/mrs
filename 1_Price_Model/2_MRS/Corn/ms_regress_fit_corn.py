"""Markov Switching Regression for Corn prices.

Translated from MS_Regress_Fit_Corn.m (based on MS_Regress toolbox by Perlin).
Uses statsmodels MarkovRegression instead of the custom MATLAB toolbox.

When run, fits the model and saves price-model outputs (k, prices, P_pr_reg)
to the inventory model case directory via the bridge module.
"""
import sys
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from bridge import extract_mrs_params, save_price_model_output

# Default output directory for inventory model case data
CASE_DIR = Path(__file__).parents[3] / '2_Inventory_Model' / '3_Case' / 'Corn'

# Price grid for discretization (matching Case_Corn.xlsx)
PRICE_GRID = np.arange(100, 850, 50, dtype=float)  # [100, 150, ..., 800]


def main():
    # Load data
    data_dir = Path(__file__).parent / 'data_Files'
    log_ret = np.loadtxt(data_dir / 'Corn_m_red.txt')

    dep = log_ret[:, 0]
    exog_raw = log_ret[:, 1:3]  # 2 explanatory variables

    # Only include exog columns that are not all-zero
    nonzero_cols = np.any(exog_raw != 0, axis=0)
    exog = exog_raw[:, nonzero_cols] if np.any(nonzero_cols) else None

    # Model specification
    # MATLAB: k=2, S=[1 0 0 1] → mean switches, variance switches
    k = 2
    switching_trend = True       # S[0]=1: intercept switches
    switching_variance = True    # S[3]=1: variance switches
    # S[1]=0, S[2]=0: exog coefficients don't switch
    switching_exog = False if exog is not None else False

    # Fit model
    mod = sm.tsa.MarkovRegression(
        dep,
        k_regimes=k,
        exog=exog,
        trend='c',
        switching_trend=switching_trend,
        switching_exog=switching_exog,
        switching_variance=switching_variance,
    )
    res = mod.fit(search_reps=20)

    # Print results
    print(res.summary())

    # Extract transition matrix (k x k)
    print('\n--- Transition Matrix ---')
    p = np.squeeze(res.regime_transition)
    if p.ndim == 1:
        # 2-regime case: params are p[0->0] and p[1->0]
        trans_mat = np.array([[p[0], 1.0 - p[1]],
                              [1.0 - p[0], p[1]]])
    else:
        trans_mat = p
    print(trans_mat)

    # Expected regime durations
    print('\n--- Expected Regime Durations ---')
    for i in range(k):
        p_ii = trans_mat[i, i]
        if p_ii < 1.0:
            duration = 1.0 / (1.0 - p_ii)
        else:
            duration = float('inf')
        print(f'  State {i + 1}: {duration:.2f} periods')

    # Plot smoothed probabilities
    fig, axes = plt.subplots(k + 1, 1, figsize=(10, 2.5 * (k + 1)), sharex=True)

    axes[0].plot(dep)
    axes[0].set_title('Dependent Variable (Corn)')
    axes[0].set_ylabel('Value')

    for i in range(k):
        axes[i + 1].plot(res.smoothed_marginal_probabilities[i])
        axes[i + 1].set_title(f'Smoothed P(State {i + 1})')
        axes[i + 1].set_ylabel('Probability')
        axes[i + 1].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(data_dir.parent / 'corn_regimes.png', dpi=150)
    plt.show()

    # Save outputs for inventory model
    regime_transition, means, sigmas, smoothed = extract_mrs_params(res)
    save_price_model_output(CASE_DIR, regime_transition, means, sigmas, PRICE_GRID)
    print(f'\nRegime means: {means}')
    print(f'Regime sigmas: {sigmas}')
    print(f'Smoothed probs (last obs): {smoothed}')


if __name__ == '__main__':
    main()
