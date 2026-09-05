import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

by_config = pd.read_csv(snakemake.input.by_config, sep="\t")


def pos_limit_sort_key(value):
    return -1 if value == "all" else float(value)


pos_limit_labels = {"all": rf"$N_{{obs}}>0$"}

pos_limits = sorted(by_config["pos_limit"].unique(), key=pos_limit_sort_key)
neg_limits = sorted(by_config["neg_limit"].unique(), key=float)


def draw_boxplot(ax, values_by_group, labels, colors):
    bp = ax.boxplot(values_by_group, positions=range(len(labels)), widths=0.5, patch_artist=True, showfliers=False)
    for box, color in zip(bp["boxes"], colors):
        box.set_facecolor(color)
        box.set_alpha(0.8)
    for median in bp["medians"]:
        median.set_color("black")

    rng = np.random.default_rng(0)
    for i, values in enumerate(values_by_group):
        x = i + rng.uniform(-0.08, 0.08, size=len(values))
        ax.scatter(x, values, color="black", s=18, zorder=3)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)


fig, axes = plt.subplots(
    1, 2, figsize=(1.8 * (len(pos_limits) + len(neg_limits)) + 2, 4.5), sharey=True
)

pos_data = [by_config.loc[by_config["pos_limit"] == p, "mean_auc"].values for p in pos_limits]
neg_data = [by_config.loc[by_config["neg_limit"] == n, "mean_auc"].values for n in neg_limits]

pos_labels = [
    pos_limit_labels[p] if p in pos_limit_labels else rf"$Q_{{2.5}} > {p}$" for p in pos_limits
]
neg_labels = [rf"$N_{{tested}} > {n}$" for n in neg_limits]

draw_boxplot(axes[0], pos_data, pos_labels, ["darkorange"] * len(pos_limits))
axes[0].set_xlabel("Positive training threshold")
axes[0].set_ylabel("ROC-AUC")

draw_boxplot(axes[1], neg_data, neg_labels, ["blue"] * len(neg_limits))
axes[1].set_xlabel("Negative training threshold")

fig.set_size_inches(3*2, 4)
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
