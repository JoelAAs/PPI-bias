import os
import re
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mean_distance_support import get_cumulative_sum

dataset_labels = {"flat": "Combined", "ms": "MS", "y2h": "Y2H"}
dataset_colors = {"Combined": "blue", "MS": "darkgreen", "Y2H": "darkorange"}

fig, ax = plt.subplots(figsize=(6, 4))

for path in snakemake.input:
    dataset = re.sub(r"^jaccard_|\.csv$", "", os.path.basename(path))
    df = pd.read_csv(path, sep="\t")

    cum_df = get_cumulative_sum(df, value_column="lower_bound_pod", cumulative_columns=["jaccard"])
    cum_df["mean_jaccard"] = cum_df["sum_jaccard"] / cum_df["non_na_pairs_jaccard"]

    label = dataset_labels.get(dataset, dataset)
    ax.plot(cum_df["value"], cum_df["mean_jaccard"], label=label, color=dataset_colors.get(label))

ax.axvline(0.15, color="grey", linestyle="--")
ax.set_xlabel(r"$Q_{2.5} > x$")
ax.set_ylabel("Mean Jaccard index, GO-BP terms")
ax.legend(title="Dataset")
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150, bbox_inches="tight")
