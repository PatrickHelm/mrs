def allocate(variables):
    ch = variables['ch']                        # holding cost
    cp = variables['cp']                        # shortage penalty cost
    d_max = variables['d_max']                  # maximum demand 
    I_max = variables['I_max']                  # maximum storage capacity
    k = variables['k']                          # regime transition probabilities
    m = variables['m']                          # number of regimes
    pi = variables['pi']                        # regime belief
    P_pr_reg = variables['P_pr_reg']            # price distributions in regimes
    t_max = variables['t_max']                  # number of periods
    prices = variables['prices']                # possible prices
    normal_mu = variables['normal']['mu']          # parameter of demand distribution
    normal_sigma = variables['normal']['sigma']    # parameter of demand distribution
    poisson_lambda = variables['poisson']['lambda']# parameter of demand distribution
    nbin_r = variables['nbin']['r']                # parameter of demand distribution
    nbin_p = variables['nbin']['p']                # parameter of demand distribution

    return ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma, poisson_lambda, nbin_r, nbin_p