import os
import urllib.request

import numpy as np
import pandas as pd


def _sanity_check(corr):
    """Fail loudly rather than silently producing a wrong answer if the upstream
    file layout has changed: square, symmetric, diagonal ~= 1, values in [-1, 1]."""
    assert corr.shape[0] == corr.shape[1], f"Corr matrix is not square {corr.shape}"
    diag = np.diag(corr)
    assert np.allclose(diag, 1.0, atol=1e-3), \
        f"Diagonal not ~1 (min={diag.min():.4f}, max={diag.max():.4f})"
    assert np.nanmin(corr) >= -1.0 - 1e-6 and np.nanmax(corr) <= 1.0 + 1e-6, \
        f"Values outside [-1, 1] (min={np.nanmin(corr):.4f}, max={np.nanmax(corr):.4f})"
    n = corr.shape[0]
    corner = corr[: min(n, 500), : min(n, 500)]
    assert np.allclose(corner, corner.T, atol=1e-3, equal_nan=True), \
        f"matrix is not symmetric (checked on a {corner.shape} corner)"


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    pkl_path = snakemake.input.pkl

    print(f"Loading primary correlation matrix (pickle): {pkl_path}", file=log, flush=True)
    corr_df = pd.read_pickle(pkl_path)
    _sanity_check(corr_df.to_numpy())
    source = f"pickle:{pkl_path} (primary, 6/2024)"

    print(f"Loaded {corr_df.shape[0]} x {corr_df.shape[1]} correlation matrix from {source}",
          file=log, flush=True)

    symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t")
    symbol_to_uniprot = symbol_map.groupby("hgnc_symbol")["uniprot_id"].apply(list).to_dict()

    genes_in_universe = [g for g in corr_df.index if g in symbol_to_uniprot]
    n_dropped = corr_df.shape[0] - len(genes_in_universe)
    print(f"{len(genes_in_universe)}/{corr_df.shape[0]} matrix genes map to a "
          f"POD-universe protein ({n_dropped} dropped, not imputed)", file=log, flush=True)

    corr = corr_df.loc[genes_in_universe, genes_in_universe].to_numpy(dtype=np.float32)
    print(f"Subset correlation matrix: {corr.shape}, {corr.nbytes / 1e6:.1f} MB",
          file=log, flush=True)

    np.save(snakemake.output.corr_npy, corr)

    index_rows = []
    for row_i, gene in enumerate(genes_in_universe):
        for uniprot_id in symbol_to_uniprot[gene]:
            index_rows.append((uniprot_id, gene, row_i))
    gene_index_df = pd.DataFrame(index_rows, columns=["uniprot_id", "hgnc_symbol", "row_index"])
    gene_index_df.to_csv(snakemake.output.gene_index, sep="\t", index=False)

    with open(snakemake.output.metadata, "w") as w:
        w.write("key\tvalue\n")
        w.write(f"matrix_source\t{source}\n")
        w.write(f"n_genes_in_matrix\t{len(genes_in_universe)}\n")
        w.write(f"n_genes_dropped_no_uniprot_map\t{n_dropped}\n")
        w.write("coexpr_method\tpearson (precomputed by ARCHS4)\n")

    print("Done.", file=log, flush=True)
    log.close()
