import numpy as np


def combvec(*arrays):
    """Create all combinations (equivalent to MATLAB combvec transposed).

    First argument varies fastest, matching MATLAB's combvec behavior.
    """
    grids = np.meshgrid(*arrays, indexing='ij')
    return np.column_stack([g.ravel(order='F') for g in grids])


def get_mins(PI_needed, prices, t, act_pi, act_price, all_mins, inv):
    """Return minimal costs from t+1 for the current state at period t.

    Direct translation of MATLAB get_mins.m.

    Parameters
    ----------
    PI_needed : list of np.ndarray
        Reachable beliefs for each period (0-indexed: PI_needed[t] is period t+1).
    prices : np.ndarray
        Vector of possible prices.
    t : int
        Current time period (0-indexed, so t=0 is period 1).
    act_pi : float
        Current belief value.
    act_price : float
        Current price.
    all_mins : np.ndarray
        All minimum costs from period t+1.
    inv : np.ndarray
        Vector of inventory levels (0 to I_max).

    Returns
    -------
    mins : np.ndarray
    """
    if t == 0:
        # In the first period, all minimum values from period t+1=2 are needed
        return all_mins
    else:
        block = len(prices)

        # Find index of act_pi in PI_needed[t]
        indices = np.where(np.isclose(PI_needed[t], act_pi))[0]
        if len(indices) == 0:
            indices = np.where(PI_needed[t] == act_pi)[0]
        index = indices[0]  # 0-based

        # Get future beliefs block from PI_needed[t+1]
        # MATLAB: ende = index * block; start = ende - (block-1)
        # Python (0-based): start = index * block, end = start + block
        start = index * block
        end_idx = start + block
        future_PI = PI_needed[t + 1][start:end_idx]

        # Find belief matching current price
        price_mask = (prices == act_price)
        PI_val = future_PI[price_mask][0]

        # Build all state combinations for period t+1
        all_comb_future = combvec(prices, inv, PI_needed[t + 1])

        # Filter by matching belief
        index_3 = np.isclose(all_comb_future[:, 2], PI_val)
        return all_mins[index_3]


def precompute_mins_indices(PI_needed, prices, t, all_mins, inv, all_comb):
    """Precompute mins slices for all states at period t.

    For use with the @njit s×y loop. Returns a flat array of all needed mins
    concatenated, plus an index array marking where each state's slice starts/ends.

    Parameters
    ----------
    PI_needed : list of np.ndarray
        Reachable beliefs.
    prices : np.ndarray
        Possible prices.
    t : int
        Current period (0-indexed).
    all_mins : np.ndarray
        Minimum costs from period t+1.
    inv : np.ndarray
        Inventory levels.
    all_comb : np.ndarray, shape (n_states, 3)
        State combinations [price, inventory, belief] for current period.

    Returns
    -------
    mins_flat : np.ndarray
        Concatenated mins for all states.
    mins_indices : np.ndarray, shape (n_states, 2)
        [start, end] indices into mins_flat for each state.
    """
    n_states = all_comb.shape[0]
    all_mins_list = []
    mins_indices = np.zeros((n_states, 2), dtype=np.int64)
    offset = 0

    for s in range(n_states):
        act_pi = all_comb[s, 2]
        act_price = all_comb[s, 0]
        m = get_mins(PI_needed, prices, t, act_pi, act_price, all_mins, inv)
        all_mins_list.append(m)
        mins_indices[s, 0] = offset
        mins_indices[s, 1] = offset + len(m)
        offset += len(m)

    mins_flat = np.concatenate(all_mins_list)
    return mins_flat, mins_indices
