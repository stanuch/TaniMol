import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

from src.config import (
    RESULTS_DIR,
    CHEMBL_VERSION,
    TARGETS,
    MORGAN_RADIUS,
    CLUSTERING_THRESHOLD,
    ACTIVITY_TYPES,
    MIN_CONFIDENCE,
)

FIGURES_DIR = RESULTS_DIR / "figures"


def _ensure_dirs():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def export_cluster_statistics(cluster_stats: dict, path: Path | None = None) -> pd.DataFrame:
    """Export per-cluster activity statistics to CSV.

    Each row is one cluster with: centroid index, size, mean/median/std/min/max pIC50.
    Sorted by cluster size descending.
    """
    _ensure_dirs()
    path = path or RESULTS_DIR / "cluster_statistics.csv"

    df = pd.DataFrame.from_dict(cluster_stats, orient="index")
    df.index.name = "centroid_idx"
    df = df.sort_values("n", ascending=False)
    df.to_csv(path)

    print(f"Exported {len(df)} cluster statistics to {path}")
    return df


def export_activity_cliffs(
    cliffs: list[dict],
    df: pd.DataFrame,
    sali_matrix: np.ndarray,
    smiles_column: str = "canonical_smiles",
    path: Path | None = None,
) -> pd.DataFrame:
    """Export activity cliff pairs to CSV with SMILES and SALI scores.

    Each row is one cliff pair with: both molecule indices, both SMILES,
    Tanimoto similarity, |ΔpIC50|, and SALI score.
    Sorted by SALI descending.
    """
    _ensure_dirs()
    path = path or RESULTS_DIR / "activity_cliffs.csv"

    if not cliffs:
        print("No activity cliffs found — skipping export.")
        return pd.DataFrame()

    cliffs_df = pd.DataFrame(cliffs)
    cliffs_df["smiles_i"] = cliffs_df["mol_i"].map(lambda x: df.iloc[x][smiles_column])
    cliffs_df["smiles_j"] = cliffs_df["mol_j"].map(lambda x: df.iloc[x][smiles_column])
    cliffs_df["sali"] = cliffs_df.apply(
        lambda row: float(sali_matrix[row["mol_i"], row["mol_j"]]), axis=1
    )
    cliffs_df = cliffs_df.sort_values("sali", ascending=False)

    cliffs_df.to_csv(path, index=False)
    print(f"Exported {len(cliffs_df)} activity cliffs to {path}")
    return cliffs_df


def export_molecule_summary(
    df: pd.DataFrame,
    clusters: dict,
    singletons: list,
    sali_matrix: np.ndarray,
    cliffs: list[dict],
    smiles_column: str = "canonical_smiles",
    pic50_column: str = "pchembl_value",
    path: Path | None = None,
) -> pd.DataFrame:
    """Export per-molecule summary to CSV.

    Each row is one molecule with: SMILES, target, pIC50, cluster ID
    (centroid index or 'singleton'), cluster size, max SALI involving
    this molecule, and number of activity cliff pairs it appears in.
    """
    _ensure_dirs()
    path = path or RESULTS_DIR / "molecule_summary.csv"

    n_mols = len(df)

    # Build cluster assignment: molecule index -> centroid index
    mol_to_cluster = {}
    mol_to_cluster_size = {}
    for centroid, members in clusters.items():
        for mol_idx in members:
            mol_to_cluster[mol_idx] = centroid
            mol_to_cluster_size[mol_idx] = len(members)

    # Max SALI per molecule (from upper triangle to avoid double counting)
    max_sali = np.zeros(n_mols)
    for i in range(n_mols):
        row_vals = sali_matrix[i, :]
        row_vals_copy = row_vals.copy()
        row_vals_copy[i] = 0  # exclude self
        max_sali[i] = row_vals_copy.max()

    # Count cliff involvement
    cliff_counts = np.zeros(n_mols, dtype=int)
    for c in cliffs:
        cliff_counts[c["mol_i"]] += 1
        cliff_counts[c["mol_j"]] += 1

    summary = pd.DataFrame(
        {
            "mol_idx": range(n_mols),
            "smiles": df[smiles_column].values,
            "target": df["target_chembl_id"].values if "target_chembl_id" in df.columns else "N/A",
            "pic50": df[pic50_column].values,
            "cluster_id": [mol_to_cluster.get(i, "singleton") for i in range(n_mols)],
            "cluster_size": [mol_to_cluster_size.get(i, 1) for i in range(n_mols)],
            "is_singleton": [i in singletons for i in range(n_mols)],
            "max_sali": max_sali,
            "n_cliff_pairs": cliff_counts,
        }
    )

    summary = summary.sort_values("pic50", ascending=False)
    summary.to_csv(path, index=False)
    print(f"Exported {len(summary)} molecule summaries to {path}")
    return summary


def export_pipeline_summary(
    df: pd.DataFrame,
    clusters: dict,
    singletons: list,
    cliffs: list[dict],
    sali_values: np.ndarray,
    rho: float,
    p_value: float,
    fingerprint_type: str = "Morgan (ECFP4)",
    path: Path | None = None,
) -> dict:
    """Export pipeline-level metadata and key results to JSON.

    Contains: dataset size, clustering results, activity cliff counts,
    SALI statistics, Spearman correlation, and pipeline parameters.
    """
    _ensure_dirs()
    path = path or RESULTS_DIR / "pipeline_summary.json"

    nonzero_sali = sali_values[sali_values > 0]

    summary = {
        "generated_at": datetime.now().isoformat(),
        "chembl_version": CHEMBL_VERSION,
        "targets": {k: v["name"] for k, v in TARGETS.items()},
        "parameters": {
            "fingerprint_type": fingerprint_type,
            "morgan_radius": MORGAN_RADIUS,
            "clustering_threshold": CLUSTERING_THRESHOLD,
            "activity_types": ACTIVITY_TYPES,
            "min_confidence_score": MIN_CONFIDENCE,
            "cliff_sim_threshold": 0.8,
            "cliff_activity_threshold": 2.0,
        },
        "dataset": {
            "n_molecules": len(df),
            "pic50_mean": round(float(df["pchembl_value"].mean()), 3),
            "pic50_median": round(float(df["pchembl_value"].median()), 3),
            "pic50_std": round(float(df["pchembl_value"].std()), 3),
            "pic50_range": [
                round(float(df["pchembl_value"].min()), 3),
                round(float(df["pchembl_value"].max()), 3),
            ],
        },
        "clustering": {
            "n_clusters": len(clusters),
            "n_singletons": len(singletons),
            "n_clustered_molecules": sum(len(c) for c in clusters.values()),
            "largest_cluster_size": max(len(c) for c in clusters.values()) if clusters else 0,
        },
        "activity_cliffs": {
            "n_cliffs": len(cliffs),
            "max_delta_pic50": round(max(c["delta_pic50"] for c in cliffs), 3) if cliffs else 0,
        },
        "sali": {
            "max": round(float(nonzero_sali.max()), 2) if len(nonzero_sali) > 0 else 0,
            "mean": round(float(nonzero_sali.mean()), 2) if len(nonzero_sali) > 0 else 0,
            "median": round(float(np.median(nonzero_sali)), 2) if len(nonzero_sali) > 0 else 0,
            "p95": round(float(np.percentile(nonzero_sali, 95)), 2) if len(nonzero_sali) > 0 else 0,
            "p99": round(float(np.percentile(nonzero_sali, 99)), 2) if len(nonzero_sali) > 0 else 0,
            "pairs_above_50": int((nonzero_sali > 50).sum()),
        },
        "spearman_correlation": {
            "rho": round(rho, 4),
            "p_value": float(p_value),
            "sar_confirmed": rho < -0.1 and p_value < 0.05,
        },
    }

    with open(path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Exported pipeline summary to {path}")
    return summary


def export_all_figures(fig_functions: dict[str, callable], dpi: int = 300) -> None:
    """Save all figures to results/figures/ as PNG.

    Parameters
    ----------
    fig_functions : dict
        Keys are filenames (without extension), values are callables
        that each create and return a matplotlib Figure.
    dpi : int
        Resolution for saved figures.
    """
    _ensure_dirs()
    import matplotlib.pyplot as plt

    for name, func in fig_functions.items():
        func()
        fig = plt.gcf()
        path = FIGURES_DIR / f"{name}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"Saved {path}")

    print(f"\nExported {len(fig_functions)} figures to {FIGURES_DIR}")
