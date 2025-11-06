import numpy as np
from numba import njit

# @njit
def get_next_prices_markov(next_price: np.ndarray,
                           observed_price: float, 
                           last_price: float, 
                           prices: np.ndarray, 
                           t: int,
                           ) -> np.ndarray:
    # Build the same price-history counter that was used to create next_price rows.
    # next_price was filled for horizon parameter t in the producer with NoP**t rows.
    # Here we want the counter length to be NoP**t as well. For t==1 this is NoP.
    if t <= 1:
        # simple case: pick row corresponding to last_price
        for i, price in enumerate(prices):
            if last_price == price:
                result = next_price[i, :]
    else:
        # Creation of all combinations of P1, P2, ..., PN (flipped cols)
        combinations = get_all_price_combinations(prices, t-1)
        mask = (combinations[:, 0] == last_price) & (combinations[:, 1] == observed_price)
        result = next_price[mask, :]
    return result

# @njit
def get_all_price_combinations(prices, t):
    """
    Python/Numba equivalent of:

        var = prices;
        for i = 1:t
            var = combvec(var, prices);
        end
        var = fliplr(var')

    prices : 1D array, length n
    t      : integer (>= 0)

    Returns
    -------
    out : 2D array, shape (n**(t+1), t+1)
          Each row is one combination, same order as MATLAB code.
    """
    num_prices = prices.size
    length = t + 1
    # total number of combinations: num_prices^(t+1)
    total = num_prices ** length

    out = np.empty((total, length), prices.dtype)

    # Use base-n representation of the row index to pick elements
    for idx in range(total):
        tmp = idx
        for c in range(length):
            digit = tmp % num_prices      # which price to take at this position
            tmp //= num_prices
            out[idx, c] = prices[digit]

    return out