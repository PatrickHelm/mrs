import numpy as np
from numba import njit

# Distribution codes
DIST_UNIFORM = 0
DIST_NORMAL = 1
DIST_POISSON = 2
DIST_NBIN = 3
DIST_MANUAL = 4


def get_demand_dist(dist_code, d_max, normal_mu=15.0, normal_sigma=3.0,
                    poisson_lambda=1.0, nbin_r=22.5, nbin_p=0.6,
                    manual_dist=None):
    """Allocates a probability distribution of the demand.

    Parameters
    ----------
    dist_code : int
        0=uniform, 1=normal, 2=poisson, 3=nbin, 4=manual
    d_max : int
        Maximum demand value.
    manual_dist : np.ndarray, optional
        If dist_code==4 and this is provided, use it directly as the demand
        distribution. Must have length d_max+1.

    Returns
    -------
    demand_dist : np.ndarray of shape (d_max+1,)
        Probability for each demand value 0, 1, ..., d_max.
    """
    real_demand = d_max + 1  # number of demand levels (0 can also be observed)
    x = np.arange(0, real_demand, dtype=np.float64)

    if dist_code == DIST_UNIFORM:
        # MATLAB: unidpdf(1:real_demand, real_demand) -> each value has prob 1/real_demand
        demand_dist = np.ones(real_demand) / real_demand

    elif dist_code == DIST_NORMAL:
        from scipy.stats import truncnorm
        a_trunc = (0 - normal_mu) / normal_sigma
        b_trunc = (d_max - normal_mu) / normal_sigma
        rv = truncnorm(a_trunc, b_trunc, loc=normal_mu, scale=normal_sigma)
        demand_dist = rv.pdf(x)
        demand_dist /= demand_dist.sum()  # ensure sums to 1

    elif dist_code == DIST_POISSON:
        from scipy.stats import poisson
        raw = poisson.pmf(x, poisson_lambda)
        demand_dist = raw / raw.sum()  # truncate & normalize

    elif dist_code == DIST_NBIN:
        from scipy.stats import nbinom
        raw = nbinom.pmf(x, nbin_r, nbin_p)
        demand_dist = raw / raw.sum()  # truncate & normalize

    elif dist_code == DIST_MANUAL:
        if manual_dist is not None:
            demand_dist = np.array(manual_dist, dtype=np.float64)
        else:
            demand_dist = np.zeros(real_demand)
            if real_demand >= 16:
                demand_dist[15] = 1.0  # deterministic demand = 15

    else:
        raise ValueError(f"Unknown dist_code: {dist_code}")

    return demand_dist.astype(np.float64)
