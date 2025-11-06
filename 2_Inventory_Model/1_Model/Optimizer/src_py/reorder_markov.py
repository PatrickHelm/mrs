import numpy as np
from numba import njit

# @njit
def reorder_markov(next_prices: np.ndarray, 
                   belief_vector: np.ndarray, 
                   maximum_inventory: int, 
                   t: int,
                   ) -> np.ndarray:
    """
    Reorder creates a large matrix with copies of 'next_prices' to fit the
    observed states in periods t.
    """
    maximum_inventory += 1  # labeling adjustment
    # Number of different prices
    num_prices = len(belief_vector)      
    # Number of inventory states for each price          
    num_inv_states = maximum_inventory * num_prices

    # determine number of times to repeat blocks
    groups = next_prices.shape[0] // num_prices
    # Each group of num_prices rows will be tiled num_inv_states times
    result = np.zeros((num_prices**t * maximum_inventory, num_prices), dtype=float)

    write_position = 0
    for g in range(groups):
        start = g * num_prices
        end = start + num_prices
        block = np.tile(next_prices[start:end, :], (num_inv_states, 1))
        end_write = write_position + block.shape[0]
        result[write_position:end_write, :] = block
        write_position = end_write

    return result
