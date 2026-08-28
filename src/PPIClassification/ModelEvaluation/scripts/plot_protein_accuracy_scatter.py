import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.protein_accuracy, sep="\t")
    df = df[df["dataset"] == snakemake.wildcards.dataset]

    train_types = ["hrni", "no"]
    train_type_labels = {"hrni": "Train: HRNI", "no": "Train: NO"}
    train_type_colors = {"hrni": "blue", "no": "darkorange"}

    fig = plt.figure(figsize=(7, 7))
    gs = fig.add_gridspec(
        2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
        wspace=0.05, hspace=0.05,
    )
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_legend.axis("off")

    for train_type in train_types:
        sub = df[df["train_type"] == train_type]
        ax.scatter(
            sub["noninteraction_ratio"],
            sub["interaction_ratio"],
            color=train_type_colors[train_type],
            s=14,
            alpha=0.6,
            edgecolor="none",
        )

    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1, zorder=0)
    ax.axvline(0.5, linestyle=":", color="grey", linewidth=1, zorder=0)
    ax.axhline(0.5, linestyle=":", color="grey", linewidth=1, zorder=0)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Ratio correct (non-interaction)")
    ax.set_ylabel("Ratio correct (interaction)")

    # Position-dodged marginal histograms: within each bin, one bar per train_type side by side.
    bins = np.linspace(0, 1, 21)
    bin_width = bins[1] - bins[0]
    sub_width = bin_width / len(train_types)
    for i, train_type in enumerate(train_types):
        sub = df[df["train_type"] == train_type]
        color = train_type_colors[train_type]

        counts_x, _ = np.histogram(sub["noninteraction_ratio"], bins=bins)
        ax_histx.bar(
            bins[:-1] + i * sub_width, counts_x, width=sub_width, align="edge",
            color=color, alpha=0.8, edgecolor="none",
        )

        counts_y, _ = np.histogram(sub["interaction_ratio"], bins=bins)
        ax_histy.barh(
            bins[:-1] + i * sub_width, counts_y, height=sub_width, align="edge",
            color=color, alpha=0.8, edgecolor="none",
        )

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax_histx.set_ylabel("Count")
    ax_histy.set_xlabel("Count")

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=train_type_colors[t], edgecolor="none", alpha=0.8, label=train_type_labels[t])
        for t in train_types
    ]
    ax_legend.legend(handles=legend_handles, loc="center", title="Training\nnegatives", frameon=False)
    dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
    dataset_title = dataset_labels.get(snakemake.wildcards.dataset, snakemake.wildcards.dataset)

    fig.suptitle(f"Per-protein prediction accuracy: {dataset_title}")

    plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
