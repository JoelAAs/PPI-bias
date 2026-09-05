import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    by_config = pd.read_csv(snakemake.input.by_config, sep="\t")

    dataset_labels = {"flat": "Combined", "ms": "MS", "y2h": "Y2H"}
    dataset_colors = {"Combined": "blue", "MS": "darkgreen", "Y2H": "darkorange"}

    datasets = sorted(by_config["dataset"].unique())

    data = [by_config.loc[by_config["dataset"] == d, "delta_auc"].values for d in datasets]

    positions = [i * 0.7 for i in range(len(datasets))]

    fig, ax = plt.subplots(figsize=(1.8 * len(datasets) + 1, 4.5))
    bp = ax.boxplot(data, positions=positions, widths=0.5, patch_artist=True, showfliers=False)
    for box, d in zip(bp["boxes"], datasets):
        box.set_facecolor(dataset_colors[dataset_labels.get(d, d)])
        box.set_edgecolor("black")
        box.set_alpha(0.8)
    for median in bp["medians"]:
        median.set_color("black")

    ax.axhline(0, linestyle="--", color="grey", linewidth=1, zorder=0)

    lo = min(v.min() for v in data)
    hi = max(v.max() for v in data)
    pad = (hi - lo) * 0.15 or 0.005
    ax.set_ylim(lo - pad, hi + pad * 2)

    ax.set_xticks(positions)
    ax.set_xticklabels([dataset_labels.get(d, d) for d in datasets])
    ax.set_ylabel("Delta ROC-AUC (train HRNI - train NO)")

    fig.set_size_inches(3, 4)
    plt.tight_layout()
    plt.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
