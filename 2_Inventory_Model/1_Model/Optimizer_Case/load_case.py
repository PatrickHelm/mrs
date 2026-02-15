"""Load case data from a directory of .txt files.

Replaces hardcoded values and Excel reads in the optimizer scripts.
"""
import numpy as np
from pathlib import Path


# --- Solution method codes ---
EXACT = 0
CEC = 1
R1 = 2
R2 = 3
OLFC = 4

METHOD_MAP = {'exact': EXACT, 'CEC': CEC, 'R1': R1, 'R2': R2, 'OLFC': OLFC}

# --- Demand distribution codes ---
DIST_MAP = {'uniform': 0, 'normal': 1, 'poisson': 2, 'nbin': 3, 'manual': 4}


def load_config(config_path):
    """Parse a key = value config file, returning a dict.

    Values are auto-typed: int if possible, then float, otherwise str.
    Lines starting with # are comments.
    """
    config = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, _, value = line.partition('=')
            key = key.strip()
            value = value.strip()
            # Try int, then float, then keep as str
            try:
                config[key] = int(value)
            except ValueError:
                try:
                    config[key] = float(value)
                except ValueError:
                    config[key] = value
    return config


def load_case_non_markov(case_dir):
    """Load all data files for the non-Markov optimizer.

    Parameters
    ----------
    case_dir : str or Path
        Directory containing config.txt, k.txt, prices.txt, P_pr_reg.txt,
        beliefs.txt, suboptimal_R1.txt, suboptimal_R2.txt, suboptimal_CEC.txt.

    Returns
    -------
    variables : dict
        Variables dict ready for backward_recursion.
    beliefs : np.ndarray
        Belief values to scan over.
    suboptimal_R1 : np.ndarray
    suboptimal_R2 : np.ndarray
    suboptimal_CEC : np.ndarray
    solution_method : str
    method_code : int
    """
    d = Path(case_dir)
    cfg = load_config(d / 'config.txt')

    k = np.loadtxt(d / 'k.txt')
    prices = np.loadtxt(d / 'prices.txt')
    P_pr_reg = np.loadtxt(d / 'P_pr_reg.txt')
    beliefs = np.loadtxt(d / 'beliefs.txt')
    if beliefs.ndim == 0:
        beliefs = beliefs.reshape(1)

    suboptimal_R1 = np.loadtxt(d / 'suboptimal_R1.txt')
    suboptimal_R2 = np.loadtxt(d / 'suboptimal_R2.txt')
    suboptimal_CEC = np.loadtxt(d / 'suboptimal_CEC.txt')
    if suboptimal_CEC.ndim == 1:
        suboptimal_CEC = suboptimal_CEC.reshape(1, -1)

    solution_method = cfg.get('solution_method', 'R1')
    method_code = METHOD_MAP[solution_method]

    # Manual demand distribution (if present)
    manual_dist_path = d / 'manual_dist.txt'
    manual_dist = (np.loadtxt(manual_dist_path)
                   if manual_dist_path.exists() else None)

    variables = {
        'ch': float(cfg['ch']),
        'cp': float(cfg['cp']),
        'd_max': int(cfg['d_max']),
        'I_max': int(cfg['I_max']),
        'k': k,
        'm': k.shape[0],
        'pi': None,  # set per iteration
        'P_pr_reg': P_pr_reg,
        't_max': int(cfg['t_max']),
        'prices': prices,
        'normal_mu': float(cfg.get('normal_mu', 15.0)),
        'normal_sigma': float(cfg.get('normal_sigma', 3.0)),
        'poisson_lambda': float(cfg.get('poisson_lambda', 1.0)),
        'nbin_r': float(cfg.get('nbin_r', 22.5)),
        'nbin_p': float(cfg.get('nbin_p', 0.6)),
    }
    if manual_dist is not None:
        variables['manual_dist'] = manual_dist

    demand_dist = cfg.get('demand_dist', 'uniform')
    dist_code = DIST_MAP[demand_dist]
    start_inventory = int(cfg.get('start_inventory', 0))

    return (variables, beliefs, suboptimal_R1, suboptimal_R2, suboptimal_CEC,
            solution_method, method_code, dist_code, start_inventory)


def load_case_markov(case_dir):
    """Load all data files for the Markov optimizer.

    Parameters
    ----------
    case_dir : str or Path
        Directory containing config.txt, k.txt, prices.txt, R1.txt, R2.txt,
        beliefs.txt, manual_dist.txt, suboptimal_R1.txt, etc.

    Returns
    -------
    variables : dict
        Variables dict ready for backward_recursion_markov.
    beliefs : np.ndarray
    suboptimal_R1 : np.ndarray
    suboptimal_R2 : np.ndarray
    suboptimal_CEC : np.ndarray
    solution_method : str
    method_code : int
    """
    d = Path(case_dir)
    cfg = load_config(d / 'config.txt')

    k = np.loadtxt(d / 'k.txt')
    prices = np.loadtxt(d / 'prices.txt')
    R1 = np.loadtxt(d / 'R1.txt')
    R2 = np.loadtxt(d / 'R2.txt')
    beliefs = np.loadtxt(d / 'beliefs.txt')
    if beliefs.ndim == 0:
        beliefs = beliefs.reshape(1)

    manual_dist_path = d / 'manual_dist.txt'
    manual_dist = (np.loadtxt(manual_dist_path)
                   if manual_dist_path.exists() else None)

    suboptimal_R1 = np.loadtxt(d / 'suboptimal_R1.txt')
    suboptimal_R2 = np.loadtxt(d / 'suboptimal_R2.txt')
    suboptimal_CEC = np.loadtxt(d / 'suboptimal_CEC.txt')
    if suboptimal_CEC.ndim == 1:
        suboptimal_CEC = suboptimal_CEC.reshape(1, -1)

    solution_method = cfg.get('solution_method', 'R2')
    method_code = METHOD_MAP[solution_method]

    variables = {
        'ch': float(cfg['ch']),
        'cp': float(cfg['cp']),
        'd_max': int(cfg['d_max']),
        'I_max': int(cfg['I_max']),
        'k': k,
        'm': k.shape[0],
        'pi': None,  # set per iteration
        'P_pr_reg': 0.0,  # not used in Markov case, but must be defined
        't_max': int(cfg['t_max']),
        'prices': prices,
        'normal_mu': float(cfg.get('normal_mu', 15.0)),
        'normal_sigma': float(cfg.get('normal_sigma', 3.0)),
        'poisson_lambda': float(cfg.get('poisson_lambda', 1.0)),
        'nbin_r': float(cfg.get('nbin_r', 22.5)),
        'nbin_p': float(cfg.get('nbin_p', 0.6)),
        'R1': R1,
        'R2': R2,
    }
    if manual_dist is not None:
        variables['manual_dist'] = manual_dist

    demand_dist = cfg.get('demand_dist', 'manual')
    dist_code = DIST_MAP[demand_dist]
    start_inventory = int(cfg.get('start_inventory', 0))

    return (variables, beliefs, suboptimal_R1, suboptimal_R2, suboptimal_CEC,
            solution_method, method_code, dist_code, start_inventory)
