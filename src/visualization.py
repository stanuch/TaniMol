import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

def plot_similarity_heatmap(matrix: np.ndarray, name: str):
    n = matrix.shape[0]

    # Clip to [0,1] to avoid floating-point negatives from float32 rounding
    dist_matrix = np.clip(1.0 - matrix, 0.0, 1.0)

    # squareform converts the (n,n) matrix to a flat condensed vector
    # Ward linkage minimises within-cluster variance
    print(f"[{name}] Computing linkage for {n:,} molecules...")
    linkage_matrix = linkage(squareform(dist_matrix), method="ward")

    # leaves_list returns molecule indices in the new order
    order = leaves_list(linkage_matrix)
    print(f"[{name}] Done. Reordering matrix...")

    # Reorder rows and columns according to the clustering
    sorted_matrix = matrix[np.ix_(order, order)]

    off_diag = matrix[np.triu_indices(n, k=1)]

    fig, ax = plt.subplots(figsize=(7, 6))
    fig.patch.set_facecolor("white")

    im = ax.imshow(
        sorted_matrix,
        cmap="viridis",
        norm=Normalize(vmin=0.0, vmax=1.0),
        aspect="auto",
        interpolation="nearest",
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Tanimoto similarity", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    cbar.outline.set_linewidth(0.5)

    ax.set_xlabel("Molecule index (clustered)", fontsize=10, labelpad=6)
    ax.set_ylabel("Molecule index (clustered)", fontsize=10, labelpad=6)
    ax.tick_params(axis="both", labelsize=9, length=3, width=0.6)

    # Tick labels still show original molecule indices at their new positions
    ticks_pos = np.linspace(0, n - 1, 6, dtype=int)
    ax.set_xticks(ticks_pos)
    ax.set_yticks(ticks_pos)
    ax.set_xticklabels(order[ticks_pos])
    ax.set_yticklabels(order[ticks_pos])

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    ax.set_title(
        f"{name} fingerprint - pairwise Tanimoto similarity  (n={n:,}, Ward clustering)",
        fontsize=11,
        fontweight="bold",
        pad=10,
    )

    stats = (
        f"mean={off_diag.mean():.3f}   "
        f"median={np.median(off_diag):.3f}   "
        f"pairs >= 0.5: {(off_diag >= 0.5).mean() * 100:.1f}%"
    )
    fig.text(
        0.5, -0.01, stats,
        ha="center", fontsize=8.5,
        color="#555555", fontstyle="italic",
        transform=fig.transFigure,
    )

    plt.tight_layout()
    plt.show()

    return order
