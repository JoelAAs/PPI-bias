"""
Module 3 - degree quadrants: 3x3 binning of deg_pos vs deg_neg (BiologicalInsights
plan §8.1-8.2). Binning is a stated interpretation of the plan's "+/-2Q" spec, not a
literal quote of it - see decision note below and plan §10 item 1.
"""
import numpy as np
import pandas as pd


def robust_z(x):
    x = np.log10(x + 1)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.zeros_like(x)
    return (x - med) / (1.4826 * mad)


def bin_z(z, cutoff):
    cat = np.full(len(z), "medium", dtype=object)
    cat[z > cutoff] = "high"
    cat[z < -cutoff] = "low"
    return cat


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    cutoff = snakemake.params.cutoff

    df = pd.read_parquet(snakemake.input.degree)
    print(f"{len(df)} proteins loaded from degree table", file=log, flush=True)
    print("## Decisions\n"
          "Binning: robust z-score on log10(degree+1), z = (x-median)/(1.4826*MAD), "
          f"bin at +/-{cutoff} (config['quadrant_bin_cutoff']). This is an "
          "interpretation of the plan's literal '+/-2Q' spec, which is undefined for "
          "non-negative degree counts on a raw scale - flagged per plan §10 item 1, "
          "not a silent choice.", file=log, flush=True)

    df["z_pos"] = robust_z(df["deg_pos"].to_numpy())
    df["z_neg"] = robust_z(df["deg_neg"].to_numpy())
    df["bin_pos"] = bin_z(df["z_pos"].to_numpy(), cutoff)
    df["bin_neg"] = bin_z(df["z_neg"].to_numpy(), cutoff)
    df["quadrant"] = df["bin_pos"] + "_pos__" + df["bin_neg"] + "_neg"

    contingency = pd.crosstab(df["bin_pos"], df["bin_neg"])
    print("3x3 contingency table (rows=bin_pos, cols=bin_neg):", file=log, flush=True)
    print(contingency.to_string(), file=log, flush=True)

    # Sanity check per plan §8.2: low/low should be large, high/high should be
    # dominated by proteins with high sum_tests.
    low_low = int(((df["bin_pos"] == "low") & (df["bin_neg"] == "low")).sum())
    high_high_mask = (df["bin_pos"] == "high") & (df["bin_neg"] == "high")
    high_high = int(high_high_mask.sum())
    median_sum_tests_hh = df.loc[high_high_mask, "sum_tests"].median() if high_high else np.nan
    median_sum_tests_all = df["sum_tests"].median()
    print(f"## Interpretation\n"
          f"Sanity check: low/low n={low_low} (expect large), high/high n={high_high} "
          f"with median sum_tests={median_sum_tests_hh} vs overall median="
          f"{median_sum_tests_all} (expect high/high >> overall). If high/high is not "
          f"dominated by high sum_tests, something is wrong with the degree table.",
          file=log, flush=True)

    df.to_csv(snakemake.output.bins, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
