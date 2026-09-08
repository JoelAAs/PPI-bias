import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

df = pd.concat(
    (pd.read_csv(f, sep="\t") for f in snakemake.input.jaccard),
    ignore_index=True,
)

dataset_labels = {"ms": "MS", "y2h": "Y2H", "flat": "Combined"}
datasets = sorted(df["dataset"].unique())  # alphabetical: flat, ms, y2h -> Combined, MS, Y2H

comparisons = [
    ("hrni_vs_hrni", "blue", "Experimental-trained vs experimental-trained"),
    ("hrni_vs_no", "darkorange", "Experimental-trained vs NO-trained"),
    ("no_vs_no", "forestgreen", "NO-trained vs NO-trained"),
]
width = 0.28
group_gap = 1.2
offsets = {name: (i - 1) * (width + 0.04) for i, (name, _, _) in enumerate(comparisons)}

fig, ax = plt.subplots(figsize=(6 + 1, 5.5))

positions, data, colors = [], [], []
for i, dataset in enumerate(datasets):
    center = i * group_gap
    subset = df[df["dataset"] == dataset]
    for name, color, _ in comparisons:
        data.append(subset.loc[subset["comparison"] == name, "jaccard"].dropna())
        positions.append(center + offsets[name])
        colors.append(color)

bp = ax.boxplot(data, positions=positions, widths=width, patch_artist=True, medianprops={"color": "black"})
for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)


def significance_label(p_value):
    if p_value < 1e-4:
        return "****"
    if p_value < 1e-3:
        return "***"
    if p_value < 1e-2:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


# HRNI-vs-NO and NO-vs-NO are unequally sized independent samples (n_permutations vs
# n_permutations-choose-2), so this is the rank-sum (Mann-Whitney U) form of the Wilcoxon test.
y_max = max((series.max() for series in data if len(series)), default=1.0)
span = y_max if y_max else 1.0
for i, dataset in enumerate(datasets):
    center = i * group_gap
    subset = df[df["dataset"] == dataset]
    hrni_no = subset.loc[subset["comparison"] == "hrni_vs_no", "jaccard"].dropna()
    no_no = subset.loc[subset["comparison"] == "no_vs_no", "jaccard"].dropna()
    if hrni_no.empty or no_no.empty:
        continue
    statistic, p_value = mannwhitneyu(hrni_no, no_no, alternative="two-sided")
    print(
        f"{dataset}\thrni_vs_no (n={len(hrni_no)}) vs no_vs_no (n={len(no_no)})\t"
        f"U={statistic:.1f}\tp={p_value:.3g}"
    )

    left = center + offsets["hrni_vs_no"]
    right = center + offsets["no_vs_no"]
    bar_y = y_max + 0.06 * span
    tick = 0.015 * span
    ax.plot([left, left, right, right], [bar_y - tick, bar_y, bar_y, bar_y - tick], color="black", lw=1.0)
    ax.text(
        (left + right) / 2,
        bar_y + 0.005 * span,
        f"{significance_label(p_value)}\np={p_value:.2g}",
        ha="center",
        va="bottom",
        fontsize=8,
    )

ax.set_ylim(top=y_max + 0.28 * span)
ax.set_xticks([i * group_gap for i in range(len(datasets))])
ax.set_xticklabels([dataset_labels.get(d, d) for d in datasets])
ax.set_ylabel("Jaccard index, prediction HRNI test set")
ax.set_xlabel("Dataset")

legend_handles = [Patch(facecolor=color, alpha=0.8, label=label) for _, color, label in comparisons]
ax.legend(handles=legend_handles, loc="best", fontsize=8)
plt.tight_layout()
plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
