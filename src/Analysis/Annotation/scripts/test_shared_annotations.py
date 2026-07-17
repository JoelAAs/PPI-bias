from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd


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


def shared(df):  # df: columns prot_a, prot_b, set_id, is_shared (0/1)
    r1 = df.loc[df.set_id == "pos", "is_shared"].mean()
    r2 = df.loc[df.set_id == "neg", "is_shared"].mean()
    return r1 - r2


def cluster_bootstrap(df, B=5000, seed=0):
    rng = np.random.default_rng(seed)
    proteins = pd.unique(df[["prot_a", "prot_b"]].values.ravel())
    idx = {p: df.index[(df.prot_a == p) | (df.prot_b == p)].to_numpy()
           for p in proteins}                      # edges touching each protein
    obs = shared(df)
    diffs = np.empty(B)
    for b in range(B):
        picked = rng.choice(proteins, size=len(proteins), replace=True)
        rows = np.concatenate([idx[p] for p in picked])   # edges follow resampled proteins
        diffs[b] = shared(df.loc[rows])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return obs, lo, hi


def _precompute_protein_stats(prot_a_idx, prot_b_idx, pos_mask, shared_matrix, n_prot):
    """
    Reduce edge-level data to per-protein sufficient statistics (computed once).

    For each protein p and annotation a, accumulate the number of edges where p
    appears (as bait or prey) and the annotation is shared, split by pos/neg set.
    This shrinks per-bootstrap-iteration work from O(n_edges) to O(n_proteins).
    """
    n_ann = shared_matrix.shape[1]
    neg_mask = ~pos_mask

    # pos is small (~6K edges): np.add.at is fine
    sm_pos = shared_matrix[pos_mask].astype(np.float64)   # (n_pos, n_ann)
    pos_sum = np.zeros((n_prot, n_ann), dtype=np.float64)
    pos_cnt = np.zeros(n_prot, dtype=np.float64)
    for role in [prot_a_idx[pos_mask], prot_b_idx[pos_mask]]:
        np.add.at(pos_sum, role, sm_pos)
        np.add.at(pos_cnt, role, 1)

    # neg is large (~19M edges): use bincount per annotation (O(n_neg) each, no huge temporaries)
    neg_bait = prot_a_idx[neg_mask]   # pre-extract once
    neg_prey = prot_b_idx[neg_mask]
    neg_sum = np.zeros((n_prot, n_ann), dtype=np.float64)
    for a in range(n_ann):
        w = shared_matrix[neg_mask, a].astype(np.float64)
        neg_sum[:, a] = (np.bincount(neg_bait, weights=w, minlength=n_prot) +
                         np.bincount(neg_prey, weights=w, minlength=n_prot))
    neg_cnt = (np.bincount(neg_bait, minlength=n_prot) +
               np.bincount(neg_prey, minlength=n_prot)).astype(np.float64)

    return pos_sum, pos_cnt, neg_sum, neg_cnt


def _log_or(r_pos, r_neg, eps=1e-10):
    """log odds ratio, clipped to avoid log(0)."""
    odds_pos = r_pos / np.maximum(1.0 - r_pos, eps)
    odds_neg = r_neg / np.maximum(1.0 - r_neg, eps)
    return np.log(np.maximum(odds_pos, eps)) - np.log(np.maximum(odds_neg, eps))


def _bootstrap_chunk(b_count, seed, n_prot, pos_sum, pos_cnt, neg_sum, neg_cnt):
    """
    Run b_count bootstrap iterations using per-protein sufficient statistics.
    O(n_prot × n_ann) per iteration instead of O(n_edges × n_ann).
    Called in a thread — arrays are shared by reference (no copying).
    """
    rng = np.random.default_rng(seed)
    n_ann = pos_sum.shape[1]
    diffs = np.empty((b_count, n_ann))
    for b in range(b_count):
        c = rng.choice(n_prot, size=n_prot, replace=True)
        bp_sum = pos_sum[c].sum(axis=0)
        bp_cnt = pos_cnt[c].sum()
        bn_sum = neg_sum[c].sum(axis=0)
        bn_cnt = neg_cnt[c].sum()
        r_pos = bp_sum / bp_cnt if bp_cnt > 0 else np.zeros(n_ann)
        r_neg = bn_sum / bn_cnt if bn_cnt > 0 else np.zeros(n_ann)
        diffs[b] = _log_or(r_pos, r_neg)
    return diffs


def cluster_bootstrap_all(df, shared_matrix, annotations, B=5000, n_workers=1, seed=0):
    """
    Cluster bootstrap over all annotations simultaneously.

    df            : DataFrame with prot_a, prot_b, set_id (reset integer index)
    shared_matrix : (n_edges, n_annotations) bool — True if both proteins carry the annotation
    annotations   : sequence of annotation labels matching shared_matrix columns

    Returns DataFrame: annotation, obs, lo, hi
    """
    prot_a = df["prot_a"].to_numpy()
    prot_b = df["prot_b"].to_numpy()
    proteins = pd.unique(np.concatenate([prot_a, prot_b]))
    prot_to_i = {p: i for i, p in enumerate(proteins)}
    n_prot = len(proteins)

    prot_a_idx = np.fromiter((prot_to_i[p] for p in prot_a), dtype=np.intp, count=len(prot_a))
    prot_b_idx = np.fromiter((prot_to_i[p] for p in prot_b), dtype=np.intp, count=len(prot_b))

    pos_mask = (df["set_id"] == "pos").to_numpy()

    rate_pos  = shared_matrix[pos_mask].mean(axis=0)
    rate_neg  = shared_matrix[~pos_mask].mean(axis=0)
    odds_pos  = rate_pos / np.maximum(1.0 - rate_pos, 1e-10)
    odds_neg  = rate_neg / np.maximum(1.0 - rate_neg, 1e-10)
    odds_ratio = odds_pos / np.maximum(odds_neg, 1e-10)

    pos_sum, pos_cnt, neg_sum, neg_cnt = _precompute_protein_stats(
        prot_a_idx, prot_b_idx, pos_mask, shared_matrix, n_prot
    )

    q, r = divmod(B, n_workers)
    chunk_sizes = [q + (1 if i < r else 0) for i in range(n_workers)]
    chunk_seeds = [seed + i * (B + 1) for i in range(n_workers)]

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [
            pool.submit(_bootstrap_chunk, b_count, s, n_prot,
                        pos_sum, pos_cnt, neg_sum, neg_cnt)
            for b_count, s in zip(chunk_sizes, chunk_seeds)
        ]
        diffs = np.vstack([f.result() for f in futures])

    # diffs contains log(OR) per bootstrap replicate; convert CI back to OR scale
    log_lo, log_hi = np.percentile(diffs, [2.5, 97.5], axis=0)
    ci_lo = np.exp(log_lo)
    ci_hi = np.exp(log_hi)

    # Two-sided p-value: fraction of bootstrap log(OR) that disagrees in sign with observed log(OR)
    log_or_obs = np.log(np.maximum(odds_ratio, 1e-10))
    p_vals = (2 * np.minimum(
        (diffs < 0).mean(axis=0),
        (diffs > 0).mean(axis=0),
    )).clip(1 / B, 1)

    # BH-FDR
    n = len(p_vals)
    order  = np.argsort(p_vals)
    ranks  = np.empty(n, dtype=int); ranks[order] = np.arange(1, n + 1)
    q_vals = np.minimum.accumulate((p_vals[order] * n / np.arange(1, n + 1))[::-1])[::-1][np.argsort(order)]

    ann_arr = np.asarray(annotations)
    return pd.DataFrame({
        "annotation":      ann_arr,
        "annotation_type": np.where(pd.Series(ann_arr).str.startswith("GO:"), "GO", "localisation"),
        "rate_pos":        rate_pos,
        "rate_neg":        rate_neg,
        "odds_pos":        odds_pos,
        "odds_neg":        odds_neg,
        "odds_ratio":      odds_ratio,
        "ci_lo":           ci_lo,
        "ci_hi":           ci_hi,
        "within_ci":       (odds_ratio >= ci_lo) & (odds_ratio <= ci_hi),
        "p_val":           p_vals,
        "q_val":           q_vals,
    })


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
