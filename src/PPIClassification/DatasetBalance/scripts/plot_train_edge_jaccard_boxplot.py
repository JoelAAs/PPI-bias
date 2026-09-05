import pandas as pd
import matplotlib.pyplot as plt

GROUPS = [
    ("hrni_negative", "Within permutation\n(HRNI negative)", "tab:blue"),
    ("no_negative", "Within permutation\n(NO negative)", "tab:green"),
    ("positive_between_config", "Between threshold config\n(positive interaction)", "darkorange"),
]


def read_jaccard(paths):
    frames = [pd.read_csv(path, sep="\t") for path in dict.fromkeys(paths)]
    return pd.concat(frames, ignore_index=True)["jaccard_index"].dropna()


data = {key: read_jaccard(getattr(snakemake.input, key)) for key, _, _ in GROUPS}

fig, ax = plt.subplots(figsize=(6, 5))

bp = ax.boxplot(
    [data[key] for key, _, _ in GROUPS],
    positions=range(len(GROUPS)),
    widths=0.55,
    patch_artist=True,
    medianprops={"color": "black"},
)
for patch, (_, _, color) in zip(bp["boxes"], GROUPS):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

ax.set_xticks(range(len(GROUPS)))
ax.set_xticklabels([f"{label}\nn={len(data[key])}" for key, label, _ in GROUPS])
ax.set_ylim(0, 1)
ax.set_ylabel("Jaccard (edge overlap)")
ax.set_title("Training-edge overlap between subsets", pad=12)

plt.tight_layout()
plt.savefig(snakemake.output.jaccard_distances, dpi=150, bbox_inches="tight")
plt.close(fig)
