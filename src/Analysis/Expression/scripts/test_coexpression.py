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
import pyarrow.dataset as ds
import pyarrow.compute as cp


#N_QUANTILE_SAMPLE_PAIRS = 200_000


def load_edges(pod_path, prot_a_col, prot_b_col, lim_col, lim_value):
    if lim_col == "lower_bound_pod":
        filter = cp.field(lim_col) >lim_value
    elif lim_col == "n_tested":
        filter = (cp.field("n_observed") == 0) & (cp.field(lim_col) >= lim_value)

    dataset = ds.dataset(pod_path)
    table = dataset.to_table(
        filter=filter,
        columns=[prot_a_col, prot_b_col]
    )
    return table.to_pandas()


def attach_coexpression(df, prot_a_col, prot_b_col, uniprot_to_row, corr):
    row_a = df[prot_a_col].map(uniprot_to_row)
    row_b = df[prot_b_col].map(uniprot_to_row)
    valid = row_a.notna() & row_b.notna()
    n_dropped = (~valid).sum()
    out = df.loc[valid].copy()
    out["coexpr"] = corr[row_a[valid].to_numpy(dtype=int), row_b[valid].to_numpy(dtype=int)]
    return out, n_dropped


def bootstrap_mean_coexp(df, n_workers):
    proteins


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    prot_a_col = snakemake.params.prot_a_col
    prot_b_col = snakemake.params.prot_b_col
    positive_max = snakemake.params.positive_max
    negative_max = snakemake.params.negative_max
    

    #rule_seed = (seed + (hash(str(snakemake.wildcards)) % 10_000)) % (2**31 - 1)

    #rng = np.random.default_rng(rule_seed)
    gene_index = pd.read_csv(snakemake.input.gene_index, sep="\t")
    uniprot_to_row = dict(zip(gene_index["uniprot_id"], gene_index["row_index"]))
    corr = np.load(snakemake.input.corr_npy)
    n_genes = corr.shape[0]

    print(f"Gene universe: {n_genes} genes ({gene_index['uniprot_id'].nunique()} uniprot IDs)",
          file=log, flush=True)


    hri_raw = load_edges(snakemake.input.pod, prot_a_col, prot_b_col, "lower_bound_pod", positive_max)
    hrni_raw = load_edges(snakemake.input.pod, prot_a_col, prot_b_col, "n_tested", negative_max)
    hri, n_drop_hri = attach_coexpression(hri_raw, prot_a_col, prot_b_col, uniprot_to_row, corr)
    hrni, n_drop_hrni = attach_coexpression(hrni_raw, prot_a_col, prot_b_col, uniprot_to_row, corr)


    hri = hri.reset_index(drop=True)
    hrni = hrni.reset_index(drop=True)
    combined = pd.concat([hri.assign(_pos=1), hrni.assign(_pos=0)], ignore_index=True)

    combined.to_csv(snakemake.output.hrni_vs_hri, sep="\t", index=None)

