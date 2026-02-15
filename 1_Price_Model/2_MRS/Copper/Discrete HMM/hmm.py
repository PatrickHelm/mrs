"""Discrete HMM training and decoding of copper prices.

Translated from HMM.m.
Uses hmmlearn CategoricalHMM instead of MATLAB's hmmtrain/hmmdecode.
"""
import numpy as np
from hmmlearn.hmm import CategoricalHMM
import matplotlib.pyplot as plt


def main():
    # 2013-2014, daily copper price buckets (8 symbols: 1-8)
    seq = [
        7, 7, 5, 5, 6, 6, 6, 6, 5, 5, 5, 5, 6, 8, 7, 7, 6, 6, 7, 7,
        7, 6, 7, 6, 7, 7, 7, 5, 6, 6, 7, 7, 7, 7, 7, 7, 8, 6, 6, 6,
        6, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 7, 7, 7, 7,
        7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8,
        7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
        6, 7, 6, 6, 6, 5, 5, 5, 5, 5, 5, 4, 3, 1, 2, 1, 1, 2, 2, 2,
        1, 1, 1, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3,
        3, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 3, 4, 4, 5, 4, 5, 5, 5, 5,
        5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4,
        4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 5, 6, 5, 5, 5, 5, 5, 5,
        4, 4, 4, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    ]
    seq = np.array(seq)

    # Convert from 1-based to 0-based symbols for hmmlearn
    seq_0based = seq - 1
    n_symbols = 8
    n_states = 2

    # Initial guesses (matching MATLAB)
    tr_guess = np.array([[0.9, 0.1],
                         [0.1, 0.9]])
    # Slightly asymmetric emission guess to break symmetry
    # (uniform guess causes EM to converge to a degenerate solution)
    rng = np.random.RandomState(42)
    emit_guess = rng.dirichlet(np.ones(n_symbols), size=n_states)

    # Fit model using Baum-Welch (EM algorithm)
    model = CategoricalHMM(
        n_components=n_states,
        n_features=n_symbols,
        n_iter=1000,
        tol=1e-6,
        init_params='',  # don't auto-initialize, use our guesses
    )
    model.startprob_ = np.array([0.5, 0.5])
    model.transmat_ = tr_guess
    model.emissionprob_ = emit_guess

    # Reshape for hmmlearn (expects 2D column vector)
    X = seq_0based.reshape(-1, 1)
    model.fit(X)

    # Results
    print('=== Estimated Transition Matrix ===')
    print(model.transmat_)

    print('\n=== Estimated Emission Matrix ===')
    print('Rows = states, Columns = symbols (1-8)')
    print(model.emissionprob_)

    # Decode: find most likely state sequence (Viterbi)
    log_prob, states = model.decode(X)
    print(f'\nLog probability of decoded sequence: {log_prob:.4f}')

    # Posterior state probabilities
    posteriors = model.predict_proba(X)

    # Plot
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    axes[0].plot(seq, 'k-', linewidth=0.5)
    axes[0].set_ylabel('Price Bucket')
    axes[0].set_title('Copper Price Buckets (1-8)')

    axes[1].plot(states + 1, 'r-', linewidth=0.5)
    axes[1].set_ylabel('State')
    axes[1].set_title('Decoded States (Viterbi)')
    axes[1].set_yticks([1, 2])

    for i in range(n_states):
        axes[2].plot(posteriors[:, i], label=f'State {i + 1}', linewidth=0.8)
    axes[2].set_ylabel('Probability')
    axes[2].set_title('Posterior State Probabilities')
    axes[2].legend()
    axes[2].set_ylim(0, 1)
    axes[2].set_xlabel('Observation')

    plt.tight_layout()
    plt.savefig('copper_discrete_hmm.png', dpi=150)
    plt.show()


if __name__ == '__main__':
    main()
