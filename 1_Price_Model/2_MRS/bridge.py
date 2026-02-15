"""Bridge between price model (MRS) and inventory model.

Extracts fitted parameters from a statsmodels MarkovRegression result
and generates the input files needed by the inventory optimizer.
"""
import numpy as np
from pathlib import Path
from scipy.stats import norm


def extract_mrs_params(res, dep_mean=0.0, dep_std=1.0):
    """Extract regime parameters from a statsmodels MarkovRegression result.

    Parameters
    ----------
    res : MarkovRegressionResults
        Fitted model result.
    dep_mean, dep_std : float
        If the dependent variable was standardized before fitting,
        provide the original mean and std to back-transform.

    Returns
    -------
    regime_transition : np.ndarray (k, k)
        Regime transition matrix with rows summing to 1.
        k[i, j] = P(next regime = j | current regime = i).
    regime_means : np.ndarray (k,)
        Back-transformed regime means.
    regime_sigmas : np.ndarray (k,)
        Back-transformed regime standard deviations.
    smoothed_probs : np.ndarray (k,)
        Smoothed regime probabilities at the last observation.
    """
    k_regimes = res.model.k_regimes
    param_names = res.model.param_names

    # Regime transition matrix
    # statsmodels: columns sum to 1 (regime_transition[i,j] = P(S_t=i|S_{t-1}=j))
    # inventory model: rows sum to 1 (k[i,j] = P(next=j|current=i))
    regime_transition = np.squeeze(res.regime_transition).T

    means = []
    sigmas = []
    for r in range(k_regimes):
        # Mean
        const_name = f'const[{r}]'
        if const_name in param_names:
            idx = param_names.index(const_name)
            means.append(res.params[idx] * dep_std + dep_mean)
        else:
            idx = param_names.index('const')
            means.append(res.params[idx] * dep_std + dep_mean)

        # Sigma
        sigma2_name = f'sigma2[{r}]'
        if sigma2_name in param_names:
            idx = param_names.index(sigma2_name)
            sigmas.append(np.sqrt(res.params[idx] * dep_std ** 2))
        else:
            idx = param_names.index('sigma2')
            sigmas.append(np.sqrt(res.params[idx] * dep_std ** 2))

    smoothed_probs = np.array([
        res.smoothed_marginal_probabilities[r][-1]
        for r in range(k_regimes)
    ])

    return regime_transition, np.array(means), np.array(sigmas), smoothed_probs


def discretize_regime_distributions(regime_means, regime_sigmas, price_grid):
    """Compute P_pr_reg by evaluating normal PDFs at grid points.

    Parameters
    ----------
    regime_means : np.ndarray (k,)
    regime_sigmas : np.ndarray (k,)
    price_grid : np.ndarray (N,)

    Returns
    -------
    P_pr_reg : np.ndarray (N, k)
        Each column is the discretized price distribution for one regime,
        normalized to sum to 1.
    """
    dx = price_grid[1] - price_grid[0]
    k = len(regime_means)
    P_pr_reg = np.zeros((len(price_grid), k))

    for r in range(k):
        probs = norm.pdf(price_grid, regime_means[r], regime_sigmas[r]) * dx
        P_pr_reg[:, r] = probs / probs.sum()

    return P_pr_reg


def save_price_model_output(output_dir, regime_transition, regime_means,
                            regime_sigmas, price_grid):
    """Save price model output files for the inventory optimizer.

    Saves k.txt, prices.txt, and P_pr_reg.txt to *output_dir*.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    np.savetxt(output_dir / 'k.txt', regime_transition)
    np.savetxt(output_dir / 'prices.txt', price_grid)

    P_pr_reg = discretize_regime_distributions(
        regime_means, regime_sigmas, price_grid)
    np.savetxt(output_dir / 'P_pr_reg.txt', P_pr_reg)

    print(f'Saved k.txt, prices.txt, P_pr_reg.txt to {output_dir}')
