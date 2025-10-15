import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform, norm, poisson, nbinom, rv_continuous

def truncate_distribution(dist, lower, upper, x):
    # Custom truncation of a distribution's PDF between lower and upper bounds
    pdf_vals = dist.pdf(x)
    pdf_vals[(x < lower) | (x > upper)] = 0
    # Normalize to ensure sum to 1 over truncated range
    total = np.sum(pdf_vals)
    if total > 0:
        pdf_vals /= total
    return pdf_vals

def get_demand_dist(dist, demand, variables):
    ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma, poisson_lambda, nbin_r, nbin_p = variables
    real_demand = demand + 1
    x = np.arange(demand + 1)

    if dist == 'uniform':
        x = np.arange(1, real_demand + 1)
        demand_dist = uniform.pdf(x, loc=1, scale=real_demand-1)
        demand_dist /= demand_dist.sum()  # Normalize to sum to 1
    elif dist == 'normal':
        pd = norm(loc=normal_mu, scale=normal_sigma)
        demand_dist = truncate_distribution(pd, 0, demand, x)
    elif dist == 'poisson':
        pd = poisson(mu=poisson_lambda)
        demand_dist = truncate_distribution(pd, 0, demand, x)
    elif dist == 'nbin':
        # nbinom in scipy uses n and p where n is number of successes, p is prob of success
        pd = nbinom(nbin_r, nbin_p)
        demand_dist = truncate_distribution(pd, 0, demand, x)
    elif dist == 'manual':
        x = np.arange(1, real_demand + 1)
        demand_dist = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    else:
        demand_dist = np.zeros_like(x)

    plt.plot(np.arange(demand + 1), demand_dist)
    plt.show()

    return demand_dist