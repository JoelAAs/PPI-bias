import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

df = pd.read_csv(snakemake.input.jaccard, sep="\t")
pivot = df.pivot_table(
    index=["dataset", "pos_limit", "neg_limit", "permutation"],
    columns="label",
    values="jaccard",
).reset_index()


def pos_limit_sort_key(value):
    return -1 if value == "all" else float(value)


markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
pos_limits = sorted(pivot["pos_limit"].unique(), key=pos_limit_sort_key)
marker_map = {pos_limit: markers[i % len(markers)] for i, pos_limit in enumerate(pos_limits)}
pos_limit_labels = {"all": "Any observation"}

datasets = sorted(pivot["dataset"].unique())
colors = plt.cm.tab10.colors
color_map = {dataset: colors[i % len(colors)] for i, dataset in enumerate(datasets)}
dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}

neg_limits = sorted(pivot["neg_limit"].unique(), key=float)

lo = min(pivot["Negative_HRNI"].min(), pivot["Interaction"].min()) - 0.02
hi = max(pivot["Negative_HRNI"].max(), pivot["Interaction"].max()) + 0.02

fig, axes = plt.subplots(
    1, len(neg_limits),
    figsize=(4 * len(neg_limits), 4),
    sharex=True, sharey=True,
    squeeze=False,
)

for col_idx, neg_limit in enumerate(neg_limits):
    ax = axes[0][col_idx]
    subset = pivot[pivot["neg_limit"] == neg_limit]
    for (dataset, pos_limit), group in subset.groupby(["dataset", "pos_limit"]):
        ax.scatter(
            group["Negative_HRNI"],
            group["Interaction"],
            color=color_map[dataset],
            marker=marker_map[pos_limit],
            edgecolor="black",
            linewidth=0.3,
            s=50,
            alpha=0.8,
        )

    ax.plot([lo, hi], [lo, hi], linestyle="--", color="grey", linewidth=1, zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(rf"$N_{{tests}} \geq {neg_limit}$")
    ax.set_xlabel("Jaccard (Negative_HRNI)")
    if col_idx == 0:
        ax.set_ylabel("Jaccard (Interaction)")

color_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[dataset], markersize=8,
           label=dataset_labels.get(dataset, dataset))
    for dataset in datasets
]
shape_handles = [
    Line2D([0], [0], marker=marker_map[pos_limit], color="black", linestyle="", markersize=8,
           label=pos_limit_labels.get(pos_limit, str(pos_limit)))
    for pos_limit in pos_limits
]
legend_dataset = fig.legend(handles=color_handles, title="Dataset", loc="lower left", bbox_to_anchor=(1.0, 0.5))
fig.add_artist(legend_dataset)
fig.legend(handles=shape_handles, title=r"$Q_{2.5}$ threshold", loc="upper left", bbox_to_anchor=(1.0, 0.5))

fig.suptitle(f"Prediction overlap (Jaccard, {snakemake.wildcards.pair_set} pairs): Interaction vs Negative_HRNI")
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
