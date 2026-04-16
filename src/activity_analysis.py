import numpy as np
from scipy.stats import spearmanr


def within_cluster_activity_distributions(clusters: dict, pic50: np.ndarray):
    results = {}
    for centroid, members in clusters.items():
        values = pic50[list(members)]
        results[centroid] = {
            "n": len(values),
            "mean": np.mean(values),
            "median": np.median(values),
            "std": np.std(values),
            "min": np.min(values),
            "max": np.max(values),
        }
    return results


def activity_cliffs(
    similarity_matrix: np.ndarray,
    delta_pic50: np.ndarray,
    sim_threshold: float = 0.8,
    activity_threshold: float = 2.0,
):
    mask = (similarity_matrix > sim_threshold) & (delta_pic50 > activity_threshold)
    np.fill_diagonal(mask, False)

    i_idx, j_idx = np.triu_indices_from(mask, k=1)
    cliff_mask = mask[i_idx, j_idx]

    cliffs = [
        {
            "mol_i": int(i),
            "mol_j": int(j),
            "similarity": float(similarity_matrix[i, j]),
            "delta_pic50": float(delta_pic50[i, j]),
        }
        for i, j in zip(i_idx[cliff_mask], j_idx[cliff_mask])
    ]

    cliffs.sort(key=lambda x: x["delta_pic50"], reverse=True)
    return cliffs


def sali(similarity_matrix: np.ndarray, delta_pic50: np.ndarray):
    distance = 1 - similarity_matrix
    sali_matrix = np.divide(
        delta_pic50, distance, out=np.zeros_like(delta_pic50), where=distance != 0
    )

    i_idx, j_idx = np.triu_indices_from(sali_matrix, k=1)
    sali_values = sali_matrix[i_idx, j_idx]

    return sali_matrix, sali_values


def similarity_activity_correlation(
    similarity_matrix: np.ndarray, delta_pic50: np.ndarray
):
    idx = np.triu_indices(similarity_matrix.shape[0], k=1)
    rho, p_value = spearmanr(similarity_matrix[idx], delta_pic50[idx])

    return float(rho), float(p_value)
