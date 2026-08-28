import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    by_config = pd.read_csv(snakemake.input.by_config, sep="\t")
    summary = pd.read_csv(snakemake.input.summary, sep="\t").set_index("dataset")

    dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
    datasets = sorted(by_config["dataset"].unique())

    data = [by_config.loc[by_config["dataset"] == d, "mean_delta_auc"].values for d in datasets]

    fig, ax = plt.subplots(figsize=(1.8 * len(datasets) + 1, 4.5))
    bp = ax.boxplot(data, positions=range(len(datasets)), widths=0.5, patch_artist=True, showfliers=False)
    for box in bp["boxes"]:
        box.set_facecolor("#7fb3d5")
        box.set_edgecolor("black")
        box.set_alpha(0.8)
    for median in bp["medians"]:
        median.set_color("black")

    rng = np.random.default_rng(0)
    for i, values in enumerate(data):
        x = i + rng.uniform(-0.08, 0.08, size=len(values))
        ax.scatter(x, values, color="black", s=18, zorder=3)

    ax.axhline(0, linestyle="--", color="grey", linewidth=1, zorder=0)

    lo = min(v.min() for v in data)
    hi = max(v.max() for v in data)
    pad = (hi - lo) * 0.15 or 0.005
    ax.set_ylim(lo - pad, hi + pad * 2)
    for i, d in enumerate(datasets):
        p = summary.loc[d, "p"]
        ax.text(i, hi + pad, f"p={p:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels([dataset_labels.get(d, d) for d in datasets])
    ax.set_ylabel("Delta ROC-AUC (train HRNI - train NO)\nper-config mean, n=6 configs")
    ax.set_title("Sister-pair delta ROC-AUC by dataset")

    plt.tight_layout()
    plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
