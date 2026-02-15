"""Markov Switching Regression for Copper prices.

Translated from MS_Regress_Fit_Copper.m (based on MS_Regress toolbox by Perlin).
Uses statsmodels MarkovRegression instead of the custom MATLAB toolbox.

Note: Data is standardized before fitting to avoid numerical issues with large
price values. Regime means and variances are back-transformed for reporting.
"""
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    # Load data
    data_dir = Path(__file__).parent / 'data_Files'
    log_ret = np.loadtxt(data_dir / 'Copper_d.txt')

    dep_raw = log_ret[:, 0]
    exog_raw = log_ret[:, 1:3]

    # Standardize to avoid numerical issues with large values
    dep_mean = dep_raw.mean()
    dep_std = dep_raw.std()
    dep = (dep_raw - dep_mean) / dep_std

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

    # Print results (standardized)
    print(res.summary())

    # Back-transform regime means and variances
    param_names = res.model.param_names
    print('\n--- Back-transformed Regime Parameters ---')
    for i in range(k):
        const_idx = param_names.index(f'const[{i}]')
        const_std = res.params[const_idx]
        const_orig = const_std * dep_std + dep_mean
        sigma2_idx = param_names.index(f'sigma2[{i}]')
        sigma2_std = res.params[sigma2_idx]
        sigma2_orig = sigma2_std * dep_std ** 2
        print(f'  Regime {i}: mean = {const_orig:.2f}, '
              f'std = {np.sqrt(sigma2_orig):.2f}')

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

    # Plot smoothed probabilities (against original data)
    fig, axes = plt.subplots(k + 1, 1, figsize=(10, 2.5 * (k + 1)), sharex=True)

    axes[0].plot(dep_raw)
    axes[0].set_title('Dependent Variable (Copper)')
    axes[0].set_ylabel('Value')

    for i in range(k):
        axes[i + 1].plot(res.smoothed_marginal_probabilities[i])
        axes[i + 1].set_title(f'Smoothed P(State {i + 1})')
        axes[i + 1].set_ylabel('Probability')
        axes[i + 1].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(data_dir.parent / 'copper_regimes.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    main()
