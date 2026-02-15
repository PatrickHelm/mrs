"""Markov Switching Regression for Zinc prices.

Translated from MS_Regress_Fit_Zinc.m (based on MS_Regress toolbox by Perlin).
Uses statsmodels MarkovRegression instead of the custom MATLAB toolbox.

When run, fits the model and saves the regime transition matrix (k.txt)
to the inventory model case directory via the bridge module.
"""
import sys
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from bridge import extract_mrs_params

# Default output directory for inventory model case data
CASE_DIR = Path(__file__).parents[3] / '2_Inventory_Model' / '3_Case' / 'Zinc'


def main():
    # Load data
    data_dir = Path(__file__).parent / 'data_Files'
    log_ret = np.loadtxt(data_dir / 'Zinc_q.txt')

    dep = log_ret[:, 0]
    exog_raw = log_ret[:, 1:3]

    # Only include exog columns that are not all-zero
    nonzero_cols = np.any(exog_raw != 0, axis=0)
    exog = exog_raw[:, nonzero_cols] if np.any(nonzero_cols) else None

    # Model specification
    # MATLAB: k=2, S=[1 0 0 1] → mean switches, variance switches
    k = 2
    switching_trend = True       # S[0]=1: intercept switches
    switching_variance = True    # S[3]=1: variance switches
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

    # Extract transition matrix
    print('\n--- Transition Matrix ---')
    trans_mat = np.squeeze(res.regime_transition)
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
    axes[0].set_title('Dependent Variable (Zinc)')
    axes[0].set_ylabel('Value')

    for i in range(k):
        axes[i + 1].plot(res.smoothed_marginal_probabilities[i])
        axes[i + 1].set_title(f'Smoothed P(State {i + 1})')
        axes[i + 1].set_ylabel('Probability')
        axes[i + 1].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(data_dir.parent / 'zinc_regimes.png', dpi=150)
    plt.show()

    # Save regime transition matrix for inventory model
    # (R1/R2 within-regime price transitions are assumed, not estimated)
    regime_transition, means, sigmas, smoothed = extract_mrs_params(res)
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    np.savetxt(CASE_DIR / 'k.txt', regime_transition)
    print(f'\nSaved k.txt to {CASE_DIR}')
    print(f'Regime means: {means}')
    print(f'Regime sigmas: {sigmas}')
    print(f'Smoothed probs (last obs): {smoothed}')


if __name__ == '__main__':
    main()
