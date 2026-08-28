import pandas as pd
import matplotlib.pyplot as plt

network_type = snakemake.wildcards.network_type


def parse_config(path):
    stem = path.split("/")[-1].removeprefix("intra_").removesuffix("_poslim-jaccard_distances.csv")
    dataset, neg_limit = stem.split(f"_{network_type}_limit_", 1)
    return dataset, neg_limit


frames = []
for path in dict.fromkeys(snakemake.input.permuted_pos):
    dataset, neg_limit = parse_config(path)
    part = pd.read_csv(path, sep="\t")
    part["dataset"] = dataset
    part["neg_limit"] = neg_limit
    frames.append(part)

df = pd.concat(frames, ignore_index=True)


def pos_limit_pair_sort_key(label):
    first = label.split("-", 1)[0]
    return -1 if first == "all" else float(first)


neg_limits = sorted(df["neg_limit"].unique(), key=float)
pos_limit_pairs = sorted(df["pos_limit"].unique(), key=pos_limit_pair_sort_key)
datasets = sorted(df["dataset"].unique())
dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
colors = plt.cm.tab10.colors
color_map = {pair: colors[i % len(colors)] for i, pair in enumerate(pos_limit_pairs)}

fig, axes = plt.subplots(
    len(neg_limits), 1,
    figsize=(4 * len(datasets), 3 * len(neg_limits)),
    sharex=True, sharey=True, squeeze=False,
)

width = 0.8 / len(pos_limit_pairs)
base_positions = list(range(len(datasets)))

for row_idx, neg_limit in enumerate(neg_limits):
    ax = axes[row_idx][0]
    subset = df[df["neg_limit"] == neg_limit]
    for pi, pos_limit_pair in enumerate(pos_limit_pairs):
        offset = (pi - (len(pos_limit_pairs) - 1) / 2) * width
        data = [
            subset.loc[
                (subset["dataset"] == d) & (subset["pos_limit"] == pos_limit_pair), "jaccard_index"
            ].dropna()
            for d in datasets
        ]
        positions = [p + offset for p in base_positions]
        bp = ax.boxplot(
            data, positions=positions, widths=width * 0.9, patch_artist=True,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color_map[pos_limit_pair])

    ax.set_xticks(base_positions)
    ax.set_xticklabels([dataset_labels.get(d, d) for d in datasets])
    ax.set_title(rf"$N_{{tests}} \geq {neg_limit}$")
    ax.set_ylabel("Jaccard (edge overlap)")

axes[-1][0].set_xlabel("Dataset")

handles = [plt.Rectangle((0, 0), 1, 1, facecolor=color_map[p]) for p in pos_limit_pairs]
fig.legend(
    handles, pos_limit_pairs,
    title="Positive threshold pair", loc="upper center", bbox_to_anchor=(0.5, 1.03),
    ncol=len(pos_limit_pairs),
)
fig.suptitle("Positive-edge overlap between adjacent positive-threshold tiers", y=1.06)
plt.tight_layout()
plt.savefig(snakemake.output.jaccard_distances, dpi=150, bbox_inches="tight")
plt.close(fig)
