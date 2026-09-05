#!/usr/bin/env python3
"""
check_log_normal.py

Tests whether the (log-)normal shape of modelled detectability is an artifact of
the GLMM's Gaussian random-effect shrinkage, or a genuine property of the data.

Logic:
  The GLMM pulls per-protein effects toward a Gaussian, which attenuates any heavy
  tail. The RAW pooled detection rate (observed / tested, per protein per role) is
  NOT pulled toward normal by any prior. So if the raw-rate distribution ALSO lacks
  a heavy tail -- especially among well-tested proteins, where shrinkage would have
  been mild anyway -- then the log-normal shape is not a shrinkage artifact.

We therefore compare, on log10 scale:
  - the raw pooled detectability distribution (all proteins, and well-tested only)
  - against a normal reference.

The row-wise input has >1e8 rows, so the per-protein counts are accumulated over
chunks instead of holding the whole frame in memory.

Input : row-wise TSV with columns  bait  prey  experiment  detection   (detection 0/1)
Output: a PNG (density overlay) and a TSV of per-protein raw rates.

Usage:
  python check_log_normal.py <row_wise_tsv> <out_png> <out_tsv> [min_tested]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

CHUNK_SIZE = 5_000_000
ROLES = ("bait", "prey")


def per_protein_rates(row_wise_path: str) -> dict:
    """Pooled observed/tested counts per protein and role, accumulated in chunks."""
    totals = {role: None for role in ROLES}
    n_rows = 0
    reader = pd.read_csv(
        row_wise_path, sep="\t",
        usecols=["bait", "prey", "detection"],
        dtype={"bait": "category", "prey": "category", "detection": "int8"},
        chunksize=CHUNK_SIZE,
    )
    for chunk in reader:
        n_rows += len(chunk)
        for role in ROLES:
            counts = chunk.groupby(role, observed=True)["detection"].agg(
                n_tested="size", n_observed="sum"
            )
            totals[role] = (
                counts if totals[role] is None
                else totals[role].add(counts, fill_value=0)
            )
        print(f"  {n_rows:,} rows", flush=True)

    print(f"Loaded {n_rows:,} rows")
    for role in ROLES:
        totals[role]["rate"] = totals[role]["n_observed"] / totals[role]["n_tested"]
        totals[role].index.name = "protein"
    return totals


def main():
    row_wise_path = sys.argv[1]
    out_png       = sys.argv[2]
    out_tsv       = sys.argv[3]
    min_tested    = int(sys.argv[4]) if len(sys.argv) > 4 else 300

    totals = per_protein_rates(row_wise_path)

    # --- per-protein raw pooled rate, each role, then averaged (matches the
    #     'average of bait & prey detectability' summary used for the degree dist) ---
    merged = pd.merge(
        totals["bait"][["rate", "n_tested"]].rename(
            columns={"rate": "bait_rate", "n_tested": "bait_n"}),
        totals["prey"][["rate", "n_tested"]].rename(
            columns={"rate": "prey_rate", "n_tested": "prey_n"}),
        left_index=True, right_index=True, how="outer",
    ).reset_index()

    # average of the two role rates where both exist, else whichever exists
    merged["avg_rate"] = merged[["bait_rate", "prey_rate"]].mean(axis=1, skipna=True)
    merged["min_n"] = merged[["bait_n", "prey_n"]].min(axis=1)  # weakest-role support

    # keep proteins with a positive rate (log scale) and some testing
    valid = merged[merged["avg_rate"] > 0].copy()
    well_tested = valid[valid["min_n"] >= min_tested].copy()

    print(f"Proteins with positive avg rate : {len(valid):,}")
    print(f"  of which well-tested (min_n>={min_tested}): {len(well_tested):,}")

    valid.to_csv(out_tsv, sep="\t", index=False)

    # --- normality of log10(raw rate), for reference (report but don't over-read at large n) ---
    def log_normality(x):
        lx = np.log10(x)
        # Shapiro caps at n=5000; subsample if larger just for the statistic
        xs = lx if len(lx) <= 5000 else np.random.default_rng(0).choice(lx, 5000, replace=False)
        W, p = stats.shapiro(xs)
        return lx, W, p

    lx_all, W_all, p_all = log_normality(valid["avg_rate"].values)
    lx_wt,  W_wt,  p_wt  = log_normality(well_tested["avg_rate"].values)

    # --- tail quantification: excess kurtosis of log10(rate) ---
    # heavy tail -> positive excess kurtosis; a clean log-normal -> ~0
    k_all = stats.kurtosis(lx_all, fisher=True)
    k_wt  = stats.kurtosis(lx_wt,  fisher=True)

    print(f"log10(avg rate)  all:        kurtosis={k_all:+.2f}  Shapiro W={W_all:.3f} p={p_all:.2e}")
    print(f"log10(avg rate)  well-tested:kurtosis={k_wt:+.2f}  Shapiro W={W_wt:.3f} p={p_wt:.2e}")

    # --- plot: density of log10(raw rate) vs fitted normal, all vs well-tested ---
    fig, ax = plt.subplots(figsize=(8, 5))
    grid = np.linspace(min(lx_all.min(), lx_wt.min()),
                       max(lx_all.max(), lx_wt.max()), 400)

    for lx, label, color in [(lx_all, f"all proteins (n={len(lx_all)})", "#888888"),
                             (lx_wt,  f"well-tested min_n>={min_tested} (n={len(lx_wt)})", "#1f77b4")]:
        kde = stats.gaussian_kde(lx)
        ax.plot(grid, kde(grid), color=color, lw=2, label=label)
        # fitted normal reference for this set
        mu, sd = lx.mean(), lx.std()
        ax.plot(grid, stats.norm.pdf(grid, mu, sd), color=color, lw=1, ls="--", alpha=0.7)

    ax.set_xlabel("log10(Average observed detectability) ")
    ax.set_ylabel("Density")
    ax.set_title("Full detectability vs log-normal reference")
    ax.legend(fontsize=8, loc="upper right")
    txt = (f"excess kurtosis (log10):\n"
           f"  all: {k_all:+.2f}\n"
           f"  well-tested: {k_wt:+.2f}")
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.9))

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"Wrote {out_png} and {out_tsv}")


if __name__ == "__main__":
    main()
