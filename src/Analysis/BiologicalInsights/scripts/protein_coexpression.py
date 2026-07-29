"""
Per-protein mean co-expression with its tested partners - the confound counterpart
to protein_degrees.py's sum_tests (total test count for that protein). Uses the SAME
edge population as sum_tests (the full POD table for this dataset/network_type, not
restricted to HRI/HRNI) so the two per-protein columns are directly comparable in the
sum_tests-vs-coexpression plot: are heavily-tested proteins more or less co-expressed
with their partners on average, simply because they were tested more?

Each edge contributes its single coexpr value to BOTH endpoints' running mean (an
edge is not "owned" by bait or prey specifically - protein_degrees.py's sum_tests
does the same double-counting via its bait+prey concat).
"""
import numpy as np
import pandas as pd
import pyarrow.dataset as ds

if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    bait_col = snakemake.params.bait_col
    prey_col = snakemake.params.prey_col

    gene_index = pd.read_csv(snakemake.input.gene_index, sep="\t")
    uniprot_to_row = dict(zip(gene_index["uniprot_id"], gene_index["row_index"]))
    corr = np.load(snakemake.input.corr_npy)

    dataset = ds.dataset(snakemake.input.pod)
    edges = dataset.to_table(columns=[bait_col, prey_col]).to_pandas()
    print(f"{len(edges)} tested edges loaded from {snakemake.input.pod}", file=log, flush=True)

    row_a = edges[bait_col].map(uniprot_to_row)
    row_b = edges[prey_col].map(uniprot_to_row)
    valid = row_a.notna() & row_b.notna()
    n_dropped = int((~valid).sum())
    edges = edges.loc[valid].copy()
    edges["coexpr"] = corr[row_a[valid].to_numpy(dtype=int), row_b[valid].to_numpy(dtype=int)]
    print(f"{len(edges)}/{n_dropped + len(edges)} edges had both endpoints in the "
          f"co-expression gene universe ({n_dropped} dropped)", file=log, flush=True)

    long = pd.concat([
        edges[[bait_col, "coexpr"]].rename(columns={bait_col: "uniprot_id"}),
        edges[[prey_col, "coexpr"]].rename(columns={prey_col: "uniprot_id"}),
    ], ignore_index=True)

    protein_coexpr = long.groupby("uniprot_id")["coexpr"].agg(
        mean_coexpr="mean", n_coexpr_partners="count"
    ).reset_index()
    print(f"{len(protein_coexpr)} proteins with >=1 co-expression-mapped tested partner",
          file=log, flush=True)

    protein_coexpr.to_csv(snakemake.output.protein_coexpr, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
