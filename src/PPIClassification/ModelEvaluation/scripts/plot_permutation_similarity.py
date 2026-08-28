import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input.similarity, sep="\t")
df["spearman"] = df[["spearman_bait", "spearman_prey"]].mean(axis=1)

label_order = ["pos", "neg_hrni", "neg_random"]
label_titles = {"pos": "Positives", "neg_hrni": "Negatives (HRNI)", "neg_random": "Negatives (random)"}


def pos_limit_sort_key(value):
    return -1 if value == "all" else float(value)


pos_limits = sorted(df["pos_limit"].unique(), key=pos_limit_sort_key)
pos_limit_labels = {"all": "Any\nobservation"}
neg_limits = sorted(df["neg_limit"].unique())
datasets = sorted(df["dataset"].unique())
dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
colors = plt.cm.tab10.colors
color_map = {dataset: colors[i % len(colors)] for i, dataset in enumerate(datasets)}


def plot_metric(value_col, ylabel, title, out_path):
    fig, axes = plt.subplots(
        len(label_order), len(neg_limits),
        figsize=(4 * len(neg_limits), 3 * len(label_order)),
        sharex=True, sharey=True, squeeze=False,
    )

    width = 0.8 / len(datasets)
    base_positions = list(range(len(pos_limits)))

    for row_idx, label in enumerate(label_order):
        for col_idx, neg_limit in enumerate(neg_limits):
            ax = axes[row_idx][col_idx]
            subset = df[(df["label"] == label) & (df["neg_limit"] == neg_limit)]
            for di, dataset in enumerate(datasets):
                offset = (di - (len(datasets) - 1) / 2) * width
                data = [
                    subset.loc[subset["pos_limit"] == pl, value_col].dropna()
                    for pl in pos_limits
                ]
                positions = [p + offset for p in base_positions]
                bp = ax.boxplot(
                    data, positions=positions, widths=width * 0.9, patch_artist=True,
                )
                for patch in bp["boxes"]:
                    patch.set_facecolor(color_map[dataset])

            ax.set_xticks(base_positions)
            ax.set_xticklabels([pos_limit_labels.get(pl, str(pl)) for pl in pos_limits])
            if row_idx == 0:
                ax.set_title(rf"$N_{{tests}} \geq {neg_limit}$")
            if col_idx == 0:
                ax.set_ylabel(f"{label_titles[label]}\n{ylabel}")
            if row_idx == len(label_order) - 1:
                ax.set_xlabel("Positive threshold")

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=color_map[d]) for d in datasets]
    fig.legend(
        handles, [dataset_labels.get(d, d) for d in datasets],
        title="Dataset", loc="upper center", bbox_to_anchor=(0.5, 1.03),
        ncol=len(datasets),
    )
    fig.suptitle(title, y=1.06)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


plot_metric(
    "jaccard_edges", "Jaccard (edge overlap)",
    "Content similarity across permutations",
    snakemake.output.jaccard_plot,
)
plot_metric(
    "spearman", "Spearman (degree)",
    "Degree similarity across permutations",
    snakemake.output.spearman_plot,
)
