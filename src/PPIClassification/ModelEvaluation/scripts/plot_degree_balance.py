import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.metrics, sep="\t")
df["delta_per_edge"] = (df["bait_degree_delta"] + df["prey_degree_delta"]) / df["num_edges"]
df["spearman"] = df[["spearman_bait", "spearman_prey"]].mean(axis=1)
df["negative_data"] = df["random"].map({True: "Non-observed", False: "HRNI"})

dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
df["dataset_label"] = df["dataset"].map(lambda d: dataset_labels.get(d, d))

datasets = sorted(df["dataset_label"].unique())
groups = ["HRNI", "Non-observed"]
colors = {"HRNI": "#e07b1a", "Non-observed": "#1a3fe0"}


def grouped_box(ax, value_col, title, ylabel):
    width = 0.35
    base_positions = list(range(len(datasets)))
    for gi, group in enumerate(groups):
        offset = (gi - 0.5) * width
        data = [
            df.loc[(df["dataset_label"] == d) & (df["negative_data"] == group), value_col].dropna()
            for d in datasets
        ]
        positions = [p + offset for p in base_positions]
        bp = ax.boxplot(data, positions=positions, widths=width * 0.9, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor(colors[group])
    ax.set_xticks(base_positions)
    ax.set_xticklabels(datasets)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Dataset")


fig, axes = plt.subplots(1, 2, figsize=(11, 5))
grouped_box(axes[0], "delta_per_edge", "Degree imbalance (magnitude)", "Degree delta per edge (bait + prey)")
grouped_box(axes[1], "spearman", "Degree agreement (rank order)", "Spearman correlation (bait/prey degree)")

handles = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[g]) for g in groups]
fig.legend(handles, groups, title="Negative data", loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2)
fig.suptitle("Train set degree balance: negatives vs. corresponding positives", y=1.12)
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
