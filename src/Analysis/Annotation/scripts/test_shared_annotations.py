from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd


def benjamini_hochberg(p_vals):
    """Benjamini-Hochberg FDR. Returns q_vals in the same order as p_vals."""
    p_vals = np.asarray(p_vals)
    n = len(p_vals)
    order = np.argsort(p_vals)
    ranks_order = np.arange(1, n + 1)
    q_vals = np.minimum.accumulate((p_vals[order] * n / ranks_order)[::-1])[::-1][np.argsort(order)]
    return np.clip(q_vals, 0, 1)


def _log_or(r_pos, r_neg, eps=1e-10):
    """log odds ratio, clipped to avoid log(0)."""
    odds_pos = r_pos / np.maximum(1.0 - r_pos, eps)
    odds_neg = r_neg / np.maximum(1.0 - r_neg, eps)
    return np.log(np.maximum(odds_pos, eps)) - np.log(np.maximum(odds_neg, eps))


def log_or_combine(pos_sum, pos_cnt, neg_sum, neg_cnt):
    """log(OR) of rate_pos vs rate_neg, for a single bootstrap replicate's sums."""
    r_pos = pos_sum / pos_cnt if pos_cnt > 0 else np.zeros_like(pos_sum)
    r_neg = neg_sum / neg_cnt if neg_cnt > 0 else np.zeros_like(neg_sum)
    return _log_or(r_pos, r_neg)


def precompute_protein_sums(prot_a_idx, prot_b_idx, pos_mask, value_matrix, n_prot):
    """
    Reduce edge-level data to per-protein sufficient statistics (computed once).

    Edges (bait-prey pairs) sharing a protein are statistically dependent, so the
    bootstrap below resamples PROTEINS with replacement (not edges); each replicate's
    statistic is derived from these per-protein sums rather than looping over edges.

    For each protein p and variable j, accumulate the summed value_matrix[:, j] over
    edges where p appears (as bait or prey), split by pos/neg edge set.

    value_matrix: (n_edges, n_vars) boolean shared-annotation indicator matrix.
    """
    n_vars = value_matrix.shape[1]
    neg_mask = ~pos_mask

    sm_pos = value_matrix[pos_mask].astype(np.float64)
    pos_sum = np.zeros((n_prot, n_vars), dtype=np.float64)
    pos_cnt = np.zeros(n_prot, dtype=np.float64)
    for role in [prot_a_idx[pos_mask], prot_b_idx[pos_mask]]:
        np.add.at(pos_sum, role, sm_pos)
        np.add.at(pos_cnt, role, 1)

    neg_bait = prot_a_idx[neg_mask]
    neg_prey = prot_b_idx[neg_mask]
    neg_sum = np.zeros((n_prot, n_vars), dtype=np.float64)
    for j in range(n_vars):
        w = value_matrix[neg_mask, j].astype(np.float64)
        neg_sum[:, j] = (np.bincount(neg_bait, weights=w, minlength=n_prot) +
                         np.bincount(neg_prey, weights=w, minlength=n_prot))
    neg_cnt = (np.bincount(neg_bait, minlength=n_prot) +
               np.bincount(neg_prey, minlength=n_prot)).astype(np.float64)

    return pos_sum, pos_cnt, neg_sum, neg_cnt


def bootstrap_chunk(b_count, seed, n_prot, pos_sum, pos_cnt, neg_sum, neg_cnt):
    """
    Run b_count bootstrap iterations using per-protein sufficient statistics.
    Called in a thread - arrays are shared by reference (no copying).
    """
    rng = np.random.default_rng(seed)
    n_vars = pos_sum.shape[1]
    diffs = np.empty((b_count, n_vars))
    for b in range(b_count):
        c = rng.choice(n_prot, size=n_prot, replace=True)
        bp_sum = pos_sum[c].sum(axis=0)
        bp_cnt = pos_cnt[c].sum()
        bn_sum = neg_sum[c].sum(axis=0)
        bn_cnt = neg_cnt[c].sum()
        diffs[b] = log_or_combine(bp_sum, bp_cnt, bn_sum, bn_cnt)
    return diffs


def cluster_bootstrap_stats(df, prot_a_col, prot_b_col, pos_col, value_matrix, labels,
                             B=5000, n_workers=1, seed=0):
    """
    Protein-level cluster bootstrap over one or more binary/rate variables
    simultaneously (effect = odds ratio of rate_pos vs rate_neg).

    df            : DataFrame with prot_a_col, prot_b_col, pos_col (boolean, "is this
                    a positive-set edge").
    value_matrix  : (n_edges, n_vars) boolean array, one column per variable in `labels`.
    labels        : sequence of variable names matching value_matrix columns.

    Returns DataFrame: label, obs_pos, obs_neg, effect, ci_lo, ci_hi, within_ci,
    p_val, q_val (BH across all variables in this call).
    """
    prot_a = df[prot_a_col].to_numpy()
    prot_b = df[prot_b_col].to_numpy()
    proteins = pd.unique(np.concatenate([prot_a, prot_b]))
    prot_to_i = {p: i for i, p in enumerate(proteins)}
    n_prot = len(proteins)

    prot_a_idx = np.fromiter((prot_to_i[p] for p in prot_a), dtype=np.intp, count=len(prot_a))
    prot_b_idx = np.fromiter((prot_to_i[p] for p in prot_b), dtype=np.intp, count=len(prot_b))
    pos_mask = df[pos_col].to_numpy().astype(bool)

    obs_pos = value_matrix[pos_mask].mean(axis=0)
    obs_neg = value_matrix[~pos_mask].mean(axis=0)
    odds_pos = obs_pos / np.maximum(1.0 - obs_pos, 1e-10)
    odds_neg = obs_neg / np.maximum(1.0 - obs_neg, 1e-10)
    obs_effect = odds_pos / np.maximum(odds_neg, 1e-10)

    pos_sum, pos_cnt, neg_sum, neg_cnt = precompute_protein_sums(
        prot_a_idx, prot_b_idx, pos_mask, value_matrix, n_prot
    )

    q, r = divmod(B, n_workers)
    chunk_sizes = [q + (1 if i < r else 0) for i in range(n_workers)]
    chunk_seeds = [seed + i * (B + 1) for i in range(n_workers)]

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [
            pool.submit(bootstrap_chunk, b_count, s, n_prot,
                        pos_sum, pos_cnt, neg_sum, neg_cnt)
            for b_count, s in zip(chunk_sizes, chunk_seeds)
        ]
        diffs = np.vstack([f.result() for f in futures])

    ci_lo, ci_hi = np.exp(np.percentile(diffs, [2.5, 97.5], axis=0))

    p_vals = (2 * np.minimum(
        (diffs < 0).mean(axis=0),
        (diffs > 0).mean(axis=0),
    )).clip(1 / B, 1)
    q_vals = benjamini_hochberg(p_vals)

    return pd.DataFrame({
        "label": np.asarray(labels),
        "obs_pos": obs_pos,
        "obs_neg": obs_neg,
        "effect": obs_effect,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "within_ci": (obs_effect >= ci_lo) & (obs_effect <= ci_hi),
        "p_val": p_vals,
        "q_val": q_vals,
    })


def cluster_bootstrap_all(df, shared_matrix, annotations, B=5000, n_workers=1, seed=0):
    """
    Reproduces the original output schema: annotation, annotation_type, rate_pos,
    rate_neg, odds_pos, odds_neg, odds_ratio, ci_lo, ci_hi, within_ci, p_val, q_val.

    df: DataFrame with columns prot_a, prot_b, set_id ("pos"/"neg").
    """
    pos_mask = (df["set_id"] == "pos")
    stats = cluster_bootstrap_stats(
        df.assign(_pos=pos_mask), "prot_a", "prot_b", "_pos", shared_matrix, annotations,
        B=B, n_workers=n_workers, seed=seed,
    )
    odds_pos = stats["obs_pos"] / np.maximum(1.0 - stats["obs_pos"], 1e-10)
    odds_neg = stats["obs_neg"] / np.maximum(1.0 - stats["obs_neg"], 1e-10)
    ann_arr = stats["label"].to_numpy()
    return pd.DataFrame({
        "annotation": ann_arr,
        "annotation_type": np.where(pd.Series(ann_arr).str.startswith("GO:"), "GO", "localisation"),
        "rate_pos": stats["obs_pos"].to_numpy(),
        "rate_neg": stats["obs_neg"].to_numpy(),
        "odds_pos": odds_pos.to_numpy(),
        "odds_neg": odds_neg.to_numpy(),
        "odds_ratio": stats["effect"].to_numpy(),
        "ci_lo": stats["ci_lo"].to_numpy(),
        "ci_hi": stats["ci_hi"].to_numpy(),
        "within_ci": stats["within_ci"].to_numpy(),
        "p_val": stats["p_val"].to_numpy(),
        "q_val": stats["q_val"].to_numpy(),
    })


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
