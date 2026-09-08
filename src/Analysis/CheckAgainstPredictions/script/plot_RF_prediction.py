# plot_RF_prediction.py
#
# FDR of RF2-PPI interaction calls as a function of the RF2-PPI score.
#
# Positives are pairs ever observed (n_observed > 0); negatives are HRNI, one set per
# n_tested threshold. At score t:
#
#   FDR(t) = FP(t) / (FP(t) + TP(t))
#
# where FP(t) / TP(t) are the negatives / positives scoring >= t. The FDR therefore
# depends on how strict the HRNI definition is, which is the point of the k sweep --
# a stricter N cutoff means fewer, better-supported negatives.
#
# Only pairs RF2-PPI actually scored are used, and the pre-selected random negatives
# (AF_subset == "NEG") are dropped so the negative set is HRNI alone.

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

cutoffs = snakemake.params.get("cutoffs", [0.5, 0.99])
n_thresholds = list(snakemake.params.n_thresholds)

# only the four columns the FDR needs: the full table is ~97M rows wide enough to
# need ~50 GB in pandas, the projection brings it under 2 GB
df = pd.read_parquet(
    snakemake.input.joined,
    columns=["n_tested", "n_observed", "RF2-PPI_score", "AF_subset"],
)
n_raw = len(df)

df = df[df["RF2-PPI_score"].notna()]
if "AF_subset" in df.columns:
    df = df[~df["AF_subset"].eq("NEG")]

tested = df["n_tested"].fillna(0)
observed = df["n_observed"].fillna(-1)

positives = np.sort(df.loc[observed > 0, "RF2-PPI_score"].to_numpy(dtype=float))
negatives = {
    k: np.sort(
        df.loc[(observed == 0) & (tested >= k), "RF2-PPI_score"].to_numpy(dtype=float)
    )
    for k in n_thresholds
}

grid = np.linspace(0, 1, 201)


def n_at_or_above(sorted_scores, thresholds):
    """Count of scores >= each threshold, via the sorted array."""
    return len(sorted_scores) - np.searchsorted(sorted_scores, thresholds, side="left")


tp = n_at_or_above(positives, grid)

# widths/alphas ramp with k so the stricter HRNI sets read as the darker lines,
# matching the HRNI styling in plot_survival
widths = np.linspace(1.2, 1.9, len(n_thresholds))
alphas = np.linspace(0.55, 1.0, len(n_thresholds))

fig, ax = plt.subplots(figsize=(7, 5))

report = []
for k, lw, alpha in zip(n_thresholds, widths, alphas):
    fp = n_at_or_above(negatives[k], grid)
    called = fp + tp
    fdr = np.divide(fp, called, out=np.full_like(grid, np.nan), where=called > 0)

    ax.plot(
        grid,
        fdr,
        color="#ff7f0e",
        ls="-",
        lw=lw,
        alpha=alpha,
        label=f"HRNI, $N\\geq{k}$ (n={len(negatives[k]):,})",
    )

    at_cutoff = {
        t: float(
            n_at_or_above(negatives[k], np.array([t]))[0]
            / max(
                n_at_or_above(negatives[k], np.array([t]))[0]
                + n_at_or_above(positives, np.array([t]))[0],
                1,
            )
        )
        for t in cutoffs
    }
    report.append((k, len(negatives[k]), at_cutoff))

for t in cutoffs:
    ax.axvline(t, color="grey", lw=0.8, ls="-.", alpha=0.7)
    ax.text(
        t,
        1.02,
        f"{t:g}",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=7,
        color="grey",
    )

# report the rate at each marked cutoff on the figure itself
lines = ["FDR at cutoff"] + [
    f"$N\\geq{k}$: " + ", ".join(f"{t:g} → {rate:.3f}" for t, rate in at_cutoff.items())
    for k, _, at_cutoff in report
]
ax.text(
    0.02,
    0.02,
    "\n".join(lines),
    transform=ax.transAxes,
    ha="left",
    va="bottom",
    fontsize=7,
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="grey", alpha=0.8),
)

ax.set_ylim(0, 1)
ax.set_xlabel("RF2-PPI interaction probability")
ax.set_ylabel("FDR: HRNI/ (HRNI + any observed)  above threshold")
ax.grid(alpha=0.25, lw=0.4)
ax.legend(loc="best", frameon=False, fontsize=8)

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"rows: {n_raw:,} -> {len(df):,} after RF2-PPI_score/AF_subset filter")
print(f"positives (n_observed > 0): {len(positives):,}")
for k, n_neg, at_cutoff in report:
    rates = ", ".join(f"FDR@{t:g}={rate:.4f}" for t, rate in at_cutoff.items())
    print(f"HRNI N>={k}: n={n_neg:,}  {rates}")
