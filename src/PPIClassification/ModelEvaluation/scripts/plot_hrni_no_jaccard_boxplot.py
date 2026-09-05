import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

df = pd.concat(
    (pd.read_csv(f, sep="\t") for f in snakemake.input.jaccard),
    ignore_index=True,
)

dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
datasets = sorted(df["dataset"].unique())  # alphabetical: flat, ms, y2h -> Combined, MS, Y2H

hrni_no_color = "darkorange"
no_no_color = "blue"
width = 0.32
group_gap = 0.8

fig, ax = plt.subplots(figsize=(2.5 * len(datasets) + 1, 5))

positions, data, colors = [], [], []
for i, dataset in enumerate(datasets):
    center = i * group_gap
    subset = df[df["dataset"] == dataset]
    data.append(subset.loc[subset["comparison"] == "hrni_vs_no", "jaccard"].dropna())
    positions.append(center - width / 2 - 0.03)
    colors.append(hrni_no_color)
    data.append(subset.loc[subset["comparison"] == "no_vs_no", "jaccard"].dropna())
    positions.append(center + width / 2 + 0.03)
    colors.append(no_no_color)

bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True, medianprops={"color": "black"})
for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

ax.set_xticks([i * group_gap for i in range(len(datasets))])
ax.set_xticklabels([dataset_labels.get(d, d) for d in datasets])
ax.set_ylabel("Jaccard index, prediction HRNI test set")
ax.set_xlabel("Dataset")

legend_handles = [
    Patch(facecolor=hrni_no_color, alpha=0.8, label="Experimental-trained vs NO-trained"),
    Patch(facecolor=no_no_color, alpha=0.8, label="NO-trained vs NO-trained"),
]
ax.legend(handles=legend_handles, loc="best")
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
