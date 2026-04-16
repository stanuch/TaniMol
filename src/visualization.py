import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from itertools import chain
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

FP_COLORS = {"Morgan": "#2196F3", "MACCS": "#4CAF50", "RDKit": "#FF5722"}
COLORS = {
    "primary": "#2A9D8F",
    "accent": "#E76F51",
    "secondary": "#F4A261",
    "neutral": "#8D99AE",
    "cliff_zone": "#FDDEDE",
    "text": "#4A4A4A",
    "annotation": "#264653",
    "annotation_border": "#B7C4CF",
}
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
        sizes = [len(v) for v in clusters[name].values()] + [1] * len(singletons[name])
        sizes = np.array(sizes)

        counts = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            counts.append(int(((sizes >= lo) & (sizes < hi)).sum()))

        color = FP_COLORS.get(name, "#888888")
        offset = (i - (len(names) - 1) / 2) * bar_width
        bars = ax.bar(
            x + offset, counts, bar_width, label=name, color=color, alpha=0.85
        )

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
    sorted_clusters = sorted(cluster_dict.values(), key=len, reverse=True)[:top_n]

    ordered_indices = list(chain.from_iterable(sorted_clusters))
    cluster_sizes = [len(c) for c in sorted_clusters]
    n_shown = len(ordered_indices)

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

    boundary = 0
    for size in cluster_sizes[:-1]:
        boundary += size
        ax.axhline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.7)
        ax.axvline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.7)

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
        [f"C{i + 1}\n(n={s})" for i, s in enumerate(cluster_sizes)],
        fontsize=7,
        rotation=45,
        ha="right",
    )
    ax.set_yticklabels(
        [f"C{i + 1} (n={s})" for i, s in enumerate(cluster_sizes)],
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
        0.5,
        -0.01,
        stats,
        ha="center",
        fontsize=STYLE["stats_fontsize"],
        color=COLORS["text"],
        fontstyle="italic",
        transform=fig.transFigure,
    )

    plt.tight_layout()
    plt.show()

    return order


def plot_cluster_activity_boxplots(
    clusters: dict,
    pic50: np.ndarray,
    top_n: int = 15,
) -> None:
    sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True)[
        :top_n
    ]

    data = [pic50[list(members)] for _, members in sorted_clusters]
    labels = [f"C{i + 1}\n(n={len(d)})" for i, d in enumerate(data)]

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    bp = ax.boxplot(
        data,
        patch_artist=True,
        labels=labels,
        widths=0.6,
        medianprops=dict(color=COLORS["accent"], linewidth=1.5),
        flierprops=dict(marker="o", markersize=3, alpha=0.5, markerfacecolor="#888"),
    )

    for patch in bp["boxes"]:
        patch.set_facecolor(COLORS["primary"])
        patch.set_alpha(0.6)

    ax.set_xlabel("Cluster", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel("pIC50", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_title(
        f"pIC50 distribution — top {top_n} clusters",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    _style_ax(ax)

    plt.tight_layout()
    plt.show()


def plot_activity_cliff_scatter(
    similarity_matrix: np.ndarray,
    delta_pic50: np.ndarray,
    sim_threshold: float = 0.8,
    activity_threshold: float = 2.0,
) -> None:
    from matplotlib.patches import Rectangle

    idx = np.triu_indices(similarity_matrix.shape[0], k=1)
    sim_vec = similarity_matrix[idx]
    delta_vec = delta_pic50[idx]

    cliff_mask = (sim_vec > sim_threshold) & (delta_vec > activity_threshold)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    # Subsample non-cliff points for readability
    non_cliff_idx = np.where(~cliff_mask)[0]
    max_points = min(20000, len(non_cliff_idx))
    rng = np.random.default_rng(42)
    sample = rng.choice(non_cliff_idx, size=max_points, replace=False)

    ax.scatter(
        sim_vec[sample],
        delta_vec[sample],
        s=2,
        alpha=0.3,
        color=COLORS["neutral"],
        rasterized=True,
        label=f"All pairs (n={len(sim_vec):,})",
    )

    # Activity cliffs on top
    if cliff_mask.sum() > 0:
        ax.scatter(
            sim_vec[cliff_mask],
            delta_vec[cliff_mask],
            s=12,
            alpha=0.9,
            color=COLORS["accent"],
            edgecolors="white",
            linewidths=0.3,
            zorder=5,
            label=f"Activity cliffs (n={cliff_mask.sum()})",
        )

    ax.axhline(
        activity_threshold,
        color=COLORS["accent"],
        linestyle="--",
        linewidth=0.8,
        alpha=0.5,
    )
    ax.axvline(
        sim_threshold, color=COLORS["accent"], linestyle="--", linewidth=0.8, alpha=0.5
    )

    # Cliff zone
    y_max = delta_vec.max() * 1.05
    rect = Rectangle(
        (sim_threshold, activity_threshold),
        1.0 - sim_threshold,
        y_max - activity_threshold,
        facecolor=COLORS["cliff_zone"],
        alpha=0.25,
        edgecolor="none",
        zorder=0,
    )
    ax.add_patch(rect)

    ax.set_xlabel("Tanimoto similarity", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel("|ΔpIC50|", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_title(
        "Activity cliff detection",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    ax.legend(fontsize=STYLE["stats_fontsize"], loc="upper left")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    _style_ax(ax)

    plt.tight_layout()
    plt.show()


def plot_sali_distribution(sali_values: np.ndarray) -> None:
    nonzero = sali_values[sali_values > 0]

    # Log-spaced bins from min to max (becuase the distribution is skewed and it looked really weird)
    bins = np.logspace(np.log10(nonzero.min()), np.log10(nonzero.max()), num=80)

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    ax.hist(
        nonzero,
        bins=bins,
        color=COLORS["primary"],
        alpha=0.7,
        edgecolor="white",
        linewidth=0.3,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    # Mark percentiles
    p95 = np.percentile(nonzero, 95)
    p99 = np.percentile(nonzero, 99)
    ax.axvline(
        p95,
        color=COLORS["secondary"],
        linestyle="--",
        linewidth=1.2,
        label=f"95th percentile = {p95:.1f}",
    )
    ax.axvline(
        p99,
        color=COLORS["accent"],
        linestyle="--",
        linewidth=1.2,
        label=f"99th percentile = {p99:.1f}",
    )

    ax.set_xlabel("SALI (log scale)", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel(
        "Number of pairs (log scale)", fontsize=STYLE["label_fontsize"], labelpad=6
    )
    ax.set_title(
        "SALI distribution (non-zero pairs)",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )
    ax.legend(fontsize=STYLE["stats_fontsize"])
    _style_ax(ax)

    stats = (
        f"max={nonzero.max():.1f}   "
        f"mean={nonzero.mean():.2f}   "
        f"median={np.median(nonzero):.2f}   "
        f"pairs > 50: {(nonzero > 50).sum()}"
    )
    fig.text(
        0.5,
        -0.01,
        stats,
        ha="center",
        fontsize=STYLE["stats_fontsize"],
        color=COLORS["text"],
        fontstyle="italic",
        transform=fig.transFigure,
    )

    plt.tight_layout()
    plt.show()


def plot_similarity_activity_density(
    similarity_matrix: np.ndarray,
    delta_pic50: np.ndarray,
    rho: float | None = None,
    p_value: float | None = None,
) -> None:
    idx = np.triu_indices(similarity_matrix.shape[0], k=1)
    sim_vec = similarity_matrix[idx]
    delta_vec = delta_pic50[idx]

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(STYLE["facecolor"])
    ax.set_facecolor(STYLE["facecolor"])

    hb = ax.hexbin(
        sim_vec,
        delta_vec,
        gridsize=60,
        cmap="YlOrRd",
        mincnt=1,
        linewidths=0.1,
    )

    cbar = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pair count", fontsize=STYLE["label_fontsize"])
    cbar.ax.tick_params(labelsize=STYLE["tick_fontsize"])
    cbar.outline.set_linewidth(0.5)

    ax.set_xlabel("Tanimoto similarity", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_ylabel("|ΔpIC50|", fontsize=STYLE["label_fontsize"], labelpad=6)
    ax.set_title(
        "Similarity–activity landscape (density)",
        fontsize=STYLE["title_fontsize"],
        fontweight="bold",
        pad=10,
    )

    if rho is not None:
        annotation = f"Spearman ρ = {rho:.4f}"
        if p_value is not None:
            annotation += f"  (p = {p_value:.2e})"
        ax.text(
            0.98,
            0.97,
            annotation,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=STYLE["stats_fontsize"],
            color=COLORS["annotation"],
            fontstyle="italic",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                alpha=0.8,
                edgecolor=COLORS["annotation_border"],
            ),
        )

    _style_ax(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)

    plt.tight_layout()
    plt.show()
