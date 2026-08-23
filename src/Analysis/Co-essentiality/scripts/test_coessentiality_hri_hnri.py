import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "Annotation" / "scripts"))
from test_shared_annotations import cluster_bootstrap_stats


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    category = pd.read_csv(snakemake.input.oe_lof_category, sep="\t", dtype={"uniprot_id": "string"})
    uniprot_to_category = dict(zip(category["uniprot_id"], category["category"]))

    hri = pd.read_csv(snakemake.input.balanced_positive, sep="\t", dtype={"bait": "string", "prey": "string"}).assign(pos=True)
    hrni = pd.read_csv(snakemake.input.balanced_negative, sep="\t", dtype={"bait": "string", "prey": "string"}).assign(pos=False)
    combined = pd.concat([hri, hrni], ignore_index=True)

    combined["bait_category"] = combined["bait"].map(uniprot_to_category)
    combined["prey_category"] = combined["prey"].map(uniprot_to_category)
    n_dropped = combined[["bait_category", "prey_category"]].isna().any(axis=1).sum()
    combined = combined.dropna(subset=["bait_category", "prey_category"]).reset_index(drop=True)
    print(f"{n_dropped} edges dropped (protein missing an oe_lof_upper category)", file=log, flush=True)
    print(f"{combined['pos'].sum()} HRI edges, {(~combined['pos']).sum()} HRNI edges", file=log, flush=True)

    both_high = ((combined["bait_category"] == "high") & (combined["prey_category"] == "high")).to_numpy()
    value_matrix = both_high.reshape(-1, 1)

    result = cluster_bootstrap_stats(
        combined, "bait", "prey", "pos", value_matrix, ["both_high"],
        B=5000, n_workers=snakemake.threads,
    ).rename(columns={"obs_pos": "rate_hri", "obs_neg": "rate_hrni", "effect": "odds_ratio"})

    print(
        f"OR both-high (HRI vs HRNI): {result['odds_ratio'].iloc[0]:.3f} "
        f"[{result['ci_lo'].iloc[0]:.3f}, {result['ci_hi'].iloc[0]:.3f}] "
        f"p={result['p_val'].iloc[0]:.4g}",
        file=log, flush=True,
    )

    result.to_csv(snakemake.output.coessentiality_OR, sep="\t", index=False)
    log.close()
