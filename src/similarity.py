"""Tanimoto similarity calculations for molecular fingerprints."""

import numpy as np
from sklearn.metrics import pairwise_distances


def calculate_similarity_matrix(fingerprints):
    """Calculate the pairwise Tanimoto similarity matrix for a set of fingerprints.

    Parameters
    ----------
    fingerprints : np.ndarray
        A 2D array of shape (n_molecules, n_bits) containing the fingerprints.

    Returns
    -------
    np.ndarray
        A symmetric 2D array of shape (n_molecules, n_molecules) containing
        the pairwise Tanimoto similarities. The diagonal should be 1.0.
    """

    # I tried using scikit-learn pairwise_distances but it was too slow,
    # also tried using just a simple loop at first but it was tragic, don't even try that.
    # Thanks to StackOverflow for the ideas

    if fingerprints.size == 0:
        return np.array([])

    fingerprints = fingerprints.astype(np.float32) # convert to float32 to save memory
    intersection = np.dot(fingerprints, fingerprints.T)
    counts = fingerprints.sum(axis=1) # sum of bits for each fingerprint
    union = counts[:, None] + counts[None, :] - intersection # union of bits for each pair of fingerprints
    
    return np.divide(intersection, union, out=np.zeros_like(intersection), where=union!=0)
            
