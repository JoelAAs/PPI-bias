import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

df_non_random = pd.read_csv(snakemake.input.non_random, sep="\t")
df_non_random["random"] = False
df_random = pd.read_csv(snakemake.input.random, sep="\t")
df_random["random"] = True
df = pd.concat([df_non_random, df_random], ignore_index=True)

def pos_limit_sort_key(value):
    return -1 if value == "all" else float(value)


markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
pos_limits = sorted(df["pos_limit"].unique(), key=pos_limit_sort_key)
marker_map = {pos_limit: markers[i % len(markers)] for i, pos_limit in enumerate(pos_limits)}
pos_limit_labels = {"all": "Any observation"}

datasets = sorted(df["dataset"].unique())
colors = plt.cm.tab10.colors
color_map = {dataset: colors[i % len(colors)] for i, dataset in enumerate(datasets)}
dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}

neg_limits = sorted(df["neg_limit"].unique(), key=float)
random_values = [False, True]
row_labels = {False: "HRNI training data", True: "Negative training data (NO)"}

lo = min(df["acc_no"].min(), df["acc_hrni"].min()) - 0.02
hi = max(df["acc_no"].max(), df["acc_hrni"].max()) + 0.02

fig, axes = plt.subplots(
    len(random_values), len(neg_limits),
    figsize=(4 * len(neg_limits), 4 * len(random_values)),
    sharex=True, sharey=True,
    squeeze=False,
)

for row_idx, random_value in enumerate(random_values):
    for col_idx, neg_limit in enumerate(neg_limits):
        ax = axes[row_idx][col_idx]
        subset = df[(df["random"] == random_value) & (df["neg_limit"] == neg_limit)]
        for (dataset, pos_limit), group in subset.groupby(["dataset", "pos_limit"]):
            ax.scatter(
                group["acc_no"],
                group["acc_hrni"],
                color=color_map[dataset],
                marker=marker_map[pos_limit],
                edgecolor="black",
                linewidth=0.3,
                s=50,
                alpha=0.8,
            )

        ax.plot([lo, hi], [lo, hi], linestyle="--", color="grey", linewidth=1, zorder=0)
        ax.axvline(0.5, linestyle=":", color="grey", linewidth=1, zorder=0)
        ax.axhline(0.5, linestyle=":", color="grey", linewidth=1, zorder=0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")

        if row_idx == 0:
            ax.set_title(rf"$N_{{tests}} \geq {neg_limit}$")
        if col_idx == 0:
            ax.set_ylabel(f"{row_labels[random_value]}\nNegative accuracy (HRNI)")
        if row_idx == len(random_values) - 1:
            ax.set_xlabel("Negative accuracy (NO)")

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

fig.suptitle("Negative prediction accuracy")
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
