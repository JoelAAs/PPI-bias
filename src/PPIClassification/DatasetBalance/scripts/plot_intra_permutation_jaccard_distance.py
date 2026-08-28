import pandas as pd
import matplotlib.pyplot as plt

network_type = snakemake.wildcards.network_type


def parse_config(path):
    stem = path.split("/")[-1].removesuffix("-jaccard_distances.csv")
    dataset, rest = stem.split(f"_{network_type}_limit_", 1)
    neg_limit, pos_limit = rest.split("_poslim_", 1)
    return dataset, neg_limit, pos_limit


frames = []
for path in dict.fromkeys(snakemake.input.permuted_pos):
    dataset, neg_limit, pos_limit = parse_config(path)
    part = pd.read_csv(path, sep="\t")
    part["dataset"] = dataset
    part["neg_limit"] = neg_limit
    part["pos_limit"] = pos_limit
    frames.append(part)

df = pd.concat(frames, ignore_index=True)


def pos_limit_sort_key(value):
    return -1 if value == "all" else float(value)


pos_limits = sorted(df["pos_limit"].unique(), key=pos_limit_sort_key)
pos_limit_labels = {"all": "Any\nobservation"}
neg_limits = sorted(df["neg_limit"].unique(), key=float)
datasets = sorted(df["dataset"].unique())
dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
colors = plt.cm.tab10.colors
color_map = {dataset: colors[i % len(colors)] for i, dataset in enumerate(datasets)}

fig, axes = plt.subplots(
    1, len(neg_limits),
    figsize=(4 * len(neg_limits), 4),
    sharex=True, sharey=True, squeeze=False,
)

width = 0.8 / len(datasets)
base_positions = list(range(len(pos_limits)))

for col_idx, neg_limit in enumerate(neg_limits):
    ax = axes[0][col_idx]
    subset = df[df["neg_limit"] == neg_limit]
    for di, dataset in enumerate(datasets):
        offset = (di - (len(datasets) - 1) / 2) * width
        data = [
            subset.loc[subset["pos_limit"] == pl, "jaccard_index"].dropna()
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
    ax.set_title(rf"$N_{{tests}} \geq {neg_limit}$")
    ax.set_xlabel("Positive threshold")
    if col_idx == 0:
        ax.set_ylabel("Jaccard (edge overlap)")

handles = [plt.Rectangle((0, 0), 1, 1, facecolor=color_map[d]) for d in datasets]
fig.legend(
    handles, [dataset_labels.get(d, d) for d in datasets],
    title="Dataset", loc="upper center", bbox_to_anchor=(0.5, 1.05),
    ncol=len(datasets),
)
fig.suptitle("Positive-edge overlap across permutations", y=1.1)
plt.tight_layout()
plt.savefig(snakemake.output.jaccard_distances, dpi=150, bbox_inches="tight")
plt.close(fig)
