import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from itertools import chain
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

FP_COLORS = {"Morgan": "#2196F3", "MACCS": "#4CAF50", "RDKit": "#FF5722"}
STYLE = {
    "facecolor": "white",
    "grid_alpha": 0.3,
    "grid_lw": 0.5,
    "spine_lw": 0.6,
    "title_fontsize": 12,
    "label_fontsize": 10,
    "tick_fontsize": 9,
    "stats_fontsize": 8.5,
}


def _style_ax(ax):
    ax.tick_params(axis="both", labelsize=STYLE["tick_fontsize"], length=3, width=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(STYLE["spine_lw"])
    ax.grid(True, alpha=STYLE["grid_alpha"], linewidth=STYLE["grid_lw"])


def plot_similarity_distribution(matrices: dict[str, np.ndarray]) -> None:
    """Overlay histogram of pairwise Tanimoto similarities for each fingerprint.

    Parameters
    ----------
    matrices : dict
        Keys are fingerprint names (e.g. 'Morgan'), values are square
        float32 numpy arrays with Tanimoto similarities on [0, 1].
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    for name, matrix in matrices.items():
        n = matrix.shape[0]
        idx = np.triu_indices(n, k=1)
        off_diag = matrix[idx]

        mean_val = off_diag.mean()
        median_val = np.median(off_diag)

        color = FP_COLORS.get(name, "#888888")
        ax.hist(
            off_diag,
            bins=100,
            density=True,
            alpha=0.55,
            color=color,
            label=(
                f"{name}  "
                f"(mean={mean_val:.3f}, "
                f"median={median_val:.3f}, "
                f"≥0.5: {(off_diag >= 0.5).mean() * 100:.1f}%)"
            ),
        )
        ax.axvline(mean_val, color=color, linewidth=1.2, linestyle="--", alpha=0.85)

    ax.set_xlabel("Tanimoto similarity", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel("Density", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_title(
        "Pairwise Tanimoto similarity distribution",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    ax.legend(fontsize=STYLE["stats_fontsize"])
    _style_ax(ax)

    plt.tight_layout()
    plt.show()


def plot_cluster_size_distribution(
    clusters: dict[str, dict],
    singletons: dict[str, list],
    bins: list[int] | None = None,
) -> None:
    """Bar chart of cluster size bins for each fingerprint, side-by-side.

    Parameters
    ----------
    clusters : dict
        Keys are fingerprint names, values are cluster_dict
        {centroid_idx: (centroid_idx, member_idx, ...)} as returned by Butina.
    singletons : dict
        Keys are fingerprint names, values are lists of singleton indices.
    bins : list of int, optional
        Right-exclusive bin edges for cluster sizes.
        Default: [1, 2, 6, 21, 51, inf]
    """
    if bins is None:
        bins = [1, 2, 6, 21, 51, int(1e9)]

    bin_labels = ["1 (singleton)", "2–5", "6–20", "21–50", ">50"]
    n_bins = len(bin_labels)
    names = list(clusters.keys())
    x = np.arange(n_bins)
    bar_width = 0.25

    fig, ax = plt.subplots(figsize=(9, 5))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    for i, name in enumerate(names):
        # Build full size list: singletons count as size 1
        sizes = [len(v) for v in clusters[name].values()] + [1] * len(singletons[name])
        sizes = np.array(sizes)

        counts = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            counts.append(int(((sizes >= lo) & (sizes < hi)).sum()))

        color = FP_COLORS.get(name, "#888888")
        offset = (i - (len(names) - 1) / 2) * bar_width
        bars = ax.bar(x + offset, counts, bar_width, label=name, color=color, alpha=0.85)

        for bar, count in zip(bars, counts):
            if count > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(counts) * 0.01,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                    color=color,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, fontsize=STYLE["tick_fontsize"])
    ax.set_xlabel("Cluster size", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel("Number of clusters", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_title(
        "Cluster size distribution by fingerprint type",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    ax.legend(fontsize=STYLE["stats_fontsize"])
    _style_ax(ax)

    plt.tight_layout()
    plt.show()


def plot_top_clusters_heatmap(
    similarity_matrix: np.ndarray,
    cluster_dict: dict,
    name: str,
    top_n: int = 10,
) -> None:
    """Heatmap of the top_n largest clusters, molecules ordered within each cluster.

    Molecules are sorted by cluster membership so each block on the diagonal
    corresponds to one cluster. Cluster boundaries are marked with white lines.

    Parameters
    ----------
    similarity_matrix : np.ndarray
        Full n×n Tanimoto similarity matrix (float32).
    cluster_dict : dict
        {centroid_idx: tuple_of_member_indices} from Butina (singletons excluded).
    name : str
        Fingerprint name used in the plot title.
    top_n : int
        Number of largest clusters to include.
    """
    # Sort clusters by descending size and take top_n
    sorted_clusters = sorted(cluster_dict.values(), key=len, reverse=True)[:top_n]

    # Build ordered index list: flatten clusters in order
    ordered_indices = list(chain.from_iterable(sorted_clusters))
    cluster_sizes = [len(c) for c in sorted_clusters]
    n_shown = len(ordered_indices)

    # Extract and reorder submatrix
    sub = similarity_matrix[np.ix_(ordered_indices, ordered_indices)]

    fig, ax = plt.subplots(figsize=(7, 6))
    fig.patch.set_facecolor(STYLE["facecolor"])

    im = ax.imshow(
        sub,
        cmap="viridis",
        norm=Normalize(vmin=0.0, vmax=1.0),
        aspect="auto",
        interpolation="nearest",
    )

    # Draw cluster boundary lines
    boundary = 0
    for size in cluster_sizes[:-1]:
        boundary += size
        ax.axhline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.7)
        ax.axvline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.7)

    # Tick at centre of each cluster block → cluster rank label
    centres = []
    pos = 0
    for size in cluster_sizes:
        centres.append(pos + size / 2)
        pos += size

    ax.set_title(
        f"{name} — top {top_n} clusters heatmap ({n_shown} molecules)",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    ax.set_xticks(centres)
    ax.set_yticks(centres)
    ax.set_xticklabels(
        [f"C{i+1}\n(n={s})" for i, s in enumerate(cluster_sizes)],
        fontsize=7,
        rotation=45,
        ha="right",
    )
    ax.set_yticklabels(
        [f"C{i+1} (n={s})" for i, s in enumerate(cluster_sizes)],
        fontsize=7,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Tanimoto similarity", fontsize=STYLE["label_fontsize"])
    cbar.ax.tick_params(labelsize=STYLE["tick_fontsize"])
    cbar.outline.set_linewidth(0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(STYLE["spine_lw"])

    plt.tight_layout()
    plt.show()


def plot_similarity_heatmap(matrix: np.ndarray, name: str) -> np.ndarray:
    """Full n×n Tanimoto heatmap sorted by UPGMA clustering.

    Use this as a quick sanity check — if the dataset has meaningful
    chemical structure, bright diagonal blocks should be visible.
    A uniform dark matrix means no clusters at the chosen threshold.

    Parameters
    ----------
    matrix : np.ndarray
        Square float32 Tanimoto similarity matrix.
    name : str
        Fingerprint name used in the plot title.

    Returns
    -------
    order : np.ndarray
        Molecule indices in UPGMA leaf order.
    """
    n = matrix.shape[0]
    dist_matrix = np.clip(1.0 - matrix, 0.0, 1.0)

    print(f"[{name}] Computing UPGMA linkage for {n:,} molecules...")
    linkage_matrix = linkage(squareform(dist_matrix), method="average")
    order = leaves_list(linkage_matrix)
    print(f"[{name}] Done. Reordering matrix...")

    sorted_matrix = matrix[np.ix_(order, order)]
    off_diag = matrix[np.triu_indices(n, k=1)]

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.patch.set_facecolor(STYLE["facecolor"])

    im = ax.imshow(
        sorted_matrix,
        cmap="viridis",
        norm=Normalize(vmin=0.0, vmax=1.0),
        aspect="auto",
        interpolation="nearest",
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Tanimoto similarity", fontsize=STYLE["label_fontsize"])
    cbar.ax.tick_params(labelsize=STYLE["tick_fontsize"])
    cbar.outline.set_linewidth(0.5)

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(STYLE["spine_lw"])

    ax.set_title(
        f"{name} — pairwise Tanimoto similarity (n={n:,}, UPGMA)",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )

    stats = (
        f"mean={off_diag.mean():.3f}   "
        f"median={np.median(off_diag):.3f}   "
        f"pairs ≥ 0.5: {(off_diag >= 0.5).mean() * 100:.1f}%"
    )
    fig.text(
        0.5, -0.01, stats,
        ha="center", fontsize=STYLE["stats_fontsize"],
        color="#555555", fontstyle="italic",
        transform=fig.transFigure,
    )

    plt.tight_layout()
    plt.show()

    return order