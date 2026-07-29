from collections import defaultdict
import numpy as np
import pandas as pd

from src.Analysis.BiologicalInsights.scripts.bootstrap import cluster_bootstrap_all


def map_uniprot_to_gene_id(uniprot_file):
    # gene_names on format: uniprot_id \t gene_id with header
    map_dict = {}
    with open(uniprot_file, 'r') as f:
        next(f)
        for line in f:
            uniprot_id, gene_id = line.strip().split('\t')
            map_dict[uniprot_id] = gene_id
    return map_dict


def build_annotation_dict(annotation_file):
    # gene_id \t annotation, one row per gene-annotation pair, no header
    annotation_dict = defaultdict(list)
    with open(annotation_file, 'r') as f:
        for line in f:
            gene_id, annotation = line.strip().split('\t')
            annotation_dict[gene_id].append(annotation)
    return annotation_dict


def build_protein_annotations(annotation_dict, uniprot_to_gene_id):
    """Map each UniProt ID to a frozenset of its annotation terms."""
    return {
        uid: frozenset(annotation_dict.get(gene_id, []))
        for uid, gene_id in uniprot_to_gene_id.items()
    }


if __name__ == "__main__":
    log = open(snakemake.log[0], 'w')

    gene_names_file  = snakemake.input.gene_names
    annotation_file  = snakemake.input.annotation
    pos_edges_file   = snakemake.input.edges_pos
    neg_edges_file   = snakemake.input.edges_neg

    bait_column = snakemake.params.bait_column
    prey_column = snakemake.params.prey_column

    print("Loading gene name mapping...", file=log, flush=True)
    uniprot_to_gene_id = map_uniprot_to_gene_id(gene_names_file)
    print(f"  {len(uniprot_to_gene_id)} UniProt → gene mappings", file=log, flush=True)

    print("Loading annotation dict...", file=log, flush=True)
    annotation_dict    = build_annotation_dict(annotation_file)
    protein_annotations = build_protein_annotations(annotation_dict, uniprot_to_gene_id)
    print(f"  {len(annotation_dict)} genes with annotations", file=log, flush=True)

    print("Loading edges...", file=log, flush=True)
    def _load_edges(path, set_id):
        df = pd.read_csv(path, sep='\t')[[bait_column, prey_column]].copy()
        # strip isoform suffixes so IDs match the canonical mapping
        df[bait_column] = df[bait_column].str.replace(r'-\d+$', '', regex=True)
        df[prey_column] = df[prey_column].str.replace(r'-\d+$', '', regex=True)
        df = df.drop_duplicates()
        df = df.rename(columns={bait_column: "prot_a", prey_column: "prot_b"})
        df["set_id"] = set_id
        return df

    pos_df   = _load_edges(pos_edges_file, "pos") # doesnt't contain any isforms, for main paper.
    neg_df   = _load_edges(neg_edges_file, "neg")
    edges_df = pd.concat([pos_df, neg_df], ignore_index=True)
    print(f"  {len(pos_df)} positive, {len(neg_df)} negative edges", file=log, flush=True)

    print("Building shared-annotation matrix...", file=log, flush=True)
    all_annotations = sorted({a for annots in protein_annotations.values() for a in annots})
    ann_to_idx      = {a: i for i, a in enumerate(all_annotations)}

    # protein × annotation boolean matrix
    proteins    = sorted(protein_annotations.keys())
    prot_to_idx = {p: i for i, p in enumerate(proteins)}
    prot_ann    = np.zeros((len(proteins), len(all_annotations)), dtype=bool)
    for p, annots in protein_annotations.items():
        pi = prot_to_idx[p]
        for a in annots:
            prot_ann[pi, ann_to_idx[a]] = True

    # per-edge: is_shared[i, j] = bait carries j AND prey carries j
    bait_pos = edges_df["prot_a"].map(prot_to_idx)
    prey_pos = edges_df["prot_b"].map(prot_to_idx)
    valid    = bait_pos.notna() & prey_pos.notna()
    bi = bait_pos.fillna(0).astype(int).to_numpy()
    pi = prey_pos.fillna(0).astype(int).to_numpy()

    shared_matrix = np.zeros((len(edges_df), len(all_annotations)), dtype=bool)
    shared_matrix[valid.to_numpy()] = (prot_ann[bi[valid]] & prot_ann[pi[valid]])
    n_missing = (~valid).sum()
    if n_missing:
        print(f"  {n_missing} edges had proteins not in annotation map (treated as no shared annotation)",
              file=log, flush=True)
    print(f"  {len(all_annotations)} annotation terms", file=log, flush=True)

    print(f"Running cluster bootstrap (B=5000, workers={snakemake.threads})...", file=log, flush=True)
    result_df = cluster_bootstrap_all(
        edges_df, shared_matrix, all_annotations,
        n_workers=snakemake.threads
    )
    print(f"  done", file=log, flush=True)

    print("Writing output...", file=log, flush=True)
    result_df.to_csv(snakemake.output[0], sep='\t', index=False)
    print("Done.", file=log, flush=True)
    log.close()
