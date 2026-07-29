"""
Module 1 - co-expression: do HRNI pairs co-express less than HRI pairs, and less than
random pairs? Three analyses (BiologicalInsights plan §6.2):
  1. categorisation: high/low/mid co-expression membership, OR of "high" HRI-vs-HRNI
  2. mean co-expression, HRNI vs HRI (raw and n_tested-matched)
  3. summed co-expression, HRNI vs size-matched random pairs (uniform null)

Decisions (per plan §10, stated rather than silently picked):
  - the random-pair null for (3) is uniform sampling only; the degree-preserving
    configuration-model null was not implemented in this pass.
  - the "large seeded random sample of gene pairs" used to set the high/low quantile
    threshold in (1) has no dedicated config key (a Monte-Carlo precision knob, not a
    scientific parameter) - fixed at N_QUANTILE_SAMPLE_PAIRS below.
"""
import numpy as np
import pandas as pd


#N_QUANTILE_SAMPLE_PAIRS = 200_000



def attach_coexpression(df, uniprot_to_row, corr):
    row_a = df["bait"].map(uniprot_to_row)
    row_b = df["prey"].map(uniprot_to_row)
    valid = row_a.notna() & row_b.notna()
    n_dropped = (~valid).sum()
    out = df.loc[valid].copy()
    out["coexpr"] = corr[row_a[valid].to_numpy(dtype=int), row_b[valid].to_numpy(dtype=int)]
    return out, n_dropped


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    gene_index = pd.read_csv(snakemake.input.gene_index, sep="\t")
    uniprot_to_row = dict(zip(gene_index["uniprot_id"], gene_index["row_index"]))
    corr = np.load(snakemake.input.corr_npy)

    hrni = pd.read_csv(snakemake.input.balanced_negative, sep="\t", dtype={"bait": "string", "prey": "string"})
    hri = pd.read_csv(snakemake.input.balanced_positive, sep="\t", dtype={"bait": "string", "prey": "string"})
    hri, n_drop_hri = attach_coexpression(hri, uniprot_to_row, corr)
    hrni, n_drop_hrni = attach_coexpression(hrni, uniprot_to_row, corr)

    hri = hri.reset_index(drop=True)
    hrni = hrni.reset_index(drop=True)
    combined = pd.concat([hri.assign(pos=1), hrni.assign(pos=0)], ignore_index=True)
    combined.to_csv(snakemake.output.hrni_vs_hri, sep="\t", index=None)

