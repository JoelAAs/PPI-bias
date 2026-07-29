"""
Module 2 - no-interaction hubs (BiologicalInsights plan §7.2-7.3): define the hub set
(high deg_neg, medium deg_pos, adequately tested) and a reference population of actual
interactors (low deg_neg, medium-or-high deg_pos, same tested-depth restriction) - hub
properties are tested for simple over/underrepresentation against this reference, not
against matched controls.

deg_neg_bin/deg_pos_bin ("high"/"medium"/"low") are precomputed once in
protein_degrees.py (shared with module 3) - see that script's docstring for the
z-score rule. This script just selects on them.
"""
import pandas as pd

if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    df = pd.read_parquet(snakemake.input.degree)
    print(f"{len(df)} proteins in the tested universe", file=log, flush=True)

    min_deg_tested = df["deg_tested"].median()
    tested_universe = df[df["deg_tested"] >= min_deg_tested].copy()
    print(f"Tested universe restricted to deg_tested >= median ({min_deg_tested}): "
          f"{len(tested_universe)} proteins", file=log, flush=True)

    hub_mask = (tested_universe["deg_neg_bin"] == "high") & (tested_universe["deg_pos_bin"] == "medium")
    hub_set = tested_universe.loc[hub_mask].reset_index(drop=True)
    print(f"Hub set: deg_neg_bin == 'high' AND deg_pos_bin == 'medium': "
          f"{len(hub_set)} proteins", file=log, flush=True)

    reference_mask = (tested_universe["deg_neg_bin"] == "low") & (tested_universe["deg_pos_bin"].isin(["medium", "high"]))
    reference_set = tested_universe.loc[reference_mask].reset_index(drop=True)
    print(f"Reference set: deg_neg_bin == 'low' AND deg_pos_bin in ('medium', 'high'): "
          f"{len(reference_set)} proteins", file=log, flush=True)

    hub_set.to_csv(snakemake.output.hub_set, sep="\t", index=False)
    reference_set.to_csv(snakemake.output.reference_set, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
