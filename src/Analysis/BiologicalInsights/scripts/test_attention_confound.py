"""
Module 3 - attention confound (BiologicalInsights plan §8.4): how much of
interaction_ratio = deg_pos / (deg_pos + deg_neg) is study attention (n_pubmed)
rather than biology? Per-protein rows are independent here (no bait/prey edge
dependency), so this uses a plain iid bootstrap, not the cluster/edge bootstrap used
elsewhere in this repo.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from src.Analysis.BiologicalInsights.scripts.bootstrap import iid_bootstrap_ci


def spearman(x, y):
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    denom = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return float((rxy - rxz * ryz) / denom) if denom > 0 else np.nan


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    seed = snakemake.params.seed
    B = snakemake.params.bootstrap_replicates

    df = pd.read_parquet(snakemake.input.degree)
    denom = df["deg_pos"] + df["deg_neg"]
    valid = denom > 0
    n_dropped = int((~valid).sum())
    df = df.loc[valid].reset_index(drop=True)
    df["interaction_ratio"] = df["deg_pos"] / denom[valid].to_numpy()
    print(f"{len(df)} proteins with deg_pos+deg_neg>0 ({n_dropped} dropped, "
          f"interaction_ratio undefined for degree-0 proteins)", file=log, flush=True)

    ratio = df["interaction_ratio"].to_numpy()
    n_pubmed = df["n_pubmed"].to_numpy()
    log_deg_tested = np.log10(df["deg_tested"].to_numpy() + 1)

    rho = spearman(ratio, n_pubmed)
    partial_rho = partial_spearman(ratio, n_pubmed, log_deg_tested)

    ci_lo, ci_hi, _ = iid_bootstrap_ci(
        len(df), lambda idx: spearman(ratio[idx], n_pubmed[idx]), B=B, seed=seed
    )
    p_ci_lo, p_ci_hi, _ = iid_bootstrap_ci(
        len(df),
        lambda idx: partial_spearman(ratio[idx], n_pubmed[idx], log_deg_tested[idx]),
        B=B, seed=seed + 1,
    )

    print(f"Spearman rho(interaction_ratio, n_pubmed) = {rho:.4f} "
          f"[{ci_lo:.4f}, {ci_hi:.4f}]", file=log, flush=True)
    print(f"Partial Spearman (controlling log10(deg_tested)) = {partial_rho:.4f} "
          f"[{p_ci_lo:.4f}, {p_ci_hi:.4f}]", file=log, flush=True)
    print("## Interpretation\n"
          f"{'High' if abs(rho) > 0.5 else 'Moderate/low'} |rho|: the quadrant "
          f"assignment is {'substantially' if abs(rho) > 0.5 else 'only partly'} a "
          "map of research attention rather than biology.", file=log, flush=True)

    result_df = pd.DataFrame([
        {"statistic": "spearman_rho", "value": rho, "ci_lo": ci_lo, "ci_hi": ci_hi, "n": len(df)},
        {"statistic": "partial_spearman_rho_controlling_log10_deg_tested",
         "value": partial_rho, "ci_lo": p_ci_lo, "ci_hi": p_ci_hi, "n": len(df)},
    ])
    result_df.to_csv(snakemake.output.attention, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
