"""Tanimoto similarity calculations for molecular fingerprints."""

import numpy as np
from pathlib import Path


def calculate_similarity_matrix(fingerprints):
    # I tried using scikit-learn pairwise_distances but it was too slow,
    # also tried using just a simple loop at first but it was tragic, don't even try that.
    # Thanks to StackOverflow for the ideas

    if fingerprints.size == 0:
        return np.array([])

    # NOTE: O(n^2) complexity. For n > 50k molecules consider
    # BulkTanimotoSimilarity or sparse matrix approach.

    fingerprints = fingerprints.astype(np.float32)  # convert to float32 to save memory
    intersection = np.dot(fingerprints, fingerprints.T)
    counts = fingerprints.sum(axis=1)  # sum of bits for each fingerprint
    union = (
        counts[:, None] + counts[None, :] - intersection
    )  # union of bits for each pair of fingerprints

    return np.divide(
        intersection, union, out=np.zeros_like(intersection), where=union != 0
    )


def save_similarity_matrix(similarity_matrix, path):
    if np.isnan(similarity_matrix).any():
        raise ValueError(
            "Warning: The similarity matrix contains NaN (Not a Number) values. Matrix is corrupted."
        )

    if similarity_matrix.shape[0] != similarity_matrix.shape[1]:
        raise ValueError(
            f"Warning: Matrix must be a square (NxN), got {similarity_matrix.shape} instead."
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, similarity_matrix)
    print(
        f"Saved {similarity_matrix.shape[0]}x{similarity_matrix.shape[1]} similarity matrix to {path}"
    )
