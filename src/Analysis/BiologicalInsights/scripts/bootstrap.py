"""
Shared protein-level cluster-bootstrap machinery.

Edges (bait-prey pairs) sharing a protein are statistically dependent, and both
positive and negative edge sets have skewed degree distributions. Every bootstrap CI
in this repo therefore resamples PROTEINS with replacement (not edges), deriving each
replicate's statistic from per-protein sufficient statistics rather than looping over
edges (O(n_prot) per iteration instead of O(n_edges)).

Refactored out of src/Analysis/Annotation/scripts/test_shared_annotations.py
(_bootstrap_chunk / _precompute_protein_stats / cluster_bootstrap_all), generalised
so a second combine function (mean difference, not just log odds ratio) can reuse the
same machinery for the co-expression module. test_shared_annotations.py now imports
cluster_bootstrap_all from here; output is numerically identical to the pre-refactor
version (verified against work_folder/analysis/shared_annotation_proportions/*.tsv).
"""
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
    """combine_fn for binary/rate variables: log(OR) of rate_pos vs rate_neg."""
    r_pos = pos_sum / pos_cnt if pos_cnt > 0 else np.zeros_like(pos_sum)
    r_neg = neg_sum / neg_cnt if neg_cnt > 0 else np.zeros_like(neg_sum)
    return _log_or(r_pos, r_neg)


def mean_diff_combine(pos_sum, pos_cnt, neg_sum, neg_cnt):
    """combine_fn for continuous variables: mean_pos - mean_neg."""
    m_pos = pos_sum / pos_cnt if pos_cnt > 0 else np.zeros_like(pos_sum)
    m_neg = neg_sum / neg_cnt if neg_cnt > 0 else np.zeros_like(neg_sum)
    return m_pos - m_neg


COMBINE_FNS = {"logor": log_or_combine, "diff": mean_diff_combine}


def precompute_protein_sums(prot_a_idx, prot_b_idx, pos_mask, value_matrix, n_prot):
    """
    Reduce edge-level data to per-protein sufficient statistics (computed once).

    For each protein p and variable j, accumulate the summed value_matrix[:, j] over
    edges where p appears (as bait or prey), split by pos/neg edge set.

    value_matrix: (n_edges, n_vars) numeric (bool or float) - e.g. shared-annotation
    indicator, or a continuous co-expression value.
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


def bootstrap_chunk(b_count, seed, n_prot, pos_sum, pos_cnt, neg_sum, neg_cnt, combine_fn):
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
        diffs[b] = combine_fn(bp_sum, bp_cnt, bn_sum, bn_cnt)
    return diffs


def cluster_bootstrap_stats(df, prot_a_col, prot_b_col, pos_col, value_matrix, labels,
                             combine="logor", B=5000, n_workers=1, seed=0):
    """
    Protein-level cluster bootstrap over one or more variables simultaneously.

    df            : DataFrame with prot_a_col, prot_b_col, pos_col (boolean, "is this
                    a positive-set edge").
    value_matrix  : (n_edges, n_vars) numeric array, one column per variable in `labels`.
    labels        : sequence of variable names matching value_matrix columns.
    combine       : "logor" (binary/rate variables, effect = odds ratio) or "diff"
                    (continuous variables, effect = mean_pos - mean_neg).

    Returns DataFrame: label, obs_pos, obs_neg, effect, ci_lo, ci_hi, within_ci,
    p_val, q_val (BH across all variables in this call).
    """
    combine_fn = COMBINE_FNS[combine]

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
    if combine == "logor":
        odds_pos = obs_pos / np.maximum(1.0 - obs_pos, 1e-10)
        odds_neg = obs_neg / np.maximum(1.0 - obs_neg, 1e-10)
        obs_effect = odds_pos / np.maximum(odds_neg, 1e-10)
    else:
        obs_effect = obs_pos - obs_neg

    pos_sum, pos_cnt, neg_sum, neg_cnt = precompute_protein_sums(
        prot_a_idx, prot_b_idx, pos_mask, value_matrix, n_prot
    )

    q, r = divmod(B, n_workers)
    chunk_sizes = [q + (1 if i < r else 0) for i in range(n_workers)]
    chunk_seeds = [seed + i * (B + 1) for i in range(n_workers)]

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [
            pool.submit(bootstrap_chunk, b_count, s, n_prot,
                        pos_sum, pos_cnt, neg_sum, neg_cnt, combine_fn)
            for b_count, s in zip(chunk_sizes, chunk_seeds)
        ]
        diffs = np.vstack([f.result() for f in futures])

    if combine == "logor":
        ci_lo, ci_hi = np.exp(np.percentile(diffs, [2.5, 97.5], axis=0))
    else:
        ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5], axis=0)

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
    Back-compat wrapper reproducing the original test_shared_annotations.py output
    schema exactly: annotation, annotation_type, rate_pos, rate_neg, odds_pos,
    odds_neg, odds_ratio, ci_lo, ci_hi, within_ci, p_val, q_val.

    df: DataFrame with columns prot_a, prot_b, set_id ("pos"/"neg").
    """
    pos_mask = (df["set_id"] == "pos")
    stats = cluster_bootstrap_stats(
        df.assign(_pos=pos_mask), "prot_a", "prot_b", "_pos", shared_matrix, annotations,
        combine="logor", B=B, n_workers=n_workers, seed=seed,
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


def precompute_protein_bin_sums(prot_a_idx, prot_b_idx, pos_mask, bin_codes, values, n_prot, n_bins):
    """
    Per-protein, per-bin sufficient statistics for a stratified mean-difference
    statistic (e.g. co-expression matched on log10(n_tested) decile). Extends
    precompute_protein_sums with a bin axis so a per-replicate weighted
    mean-diff-across-bins reduces to numpy array sums over resampled proteins
    (no pandas, no qcut, no per-replicate DataFrame copy).
    """
    neg_mask = ~pos_mask
    flat_size = n_prot * n_bins

    def accumulate(mask):
        sums = np.zeros(flat_size, dtype=np.float64)
        cnts = np.zeros(flat_size, dtype=np.float64)
        for role_idx in (prot_a_idx, prot_b_idx):
            flat = role_idx[mask] * n_bins + bin_codes[mask]
            sums += np.bincount(flat, weights=values[mask], minlength=flat_size)
            cnts += np.bincount(flat, minlength=flat_size)
        return sums.reshape(n_prot, n_bins), cnts.reshape(n_prot, n_bins)

    pos_sum, pos_cnt = accumulate(pos_mask)
    neg_sum, neg_cnt = accumulate(neg_mask)
    return pos_sum, pos_cnt, neg_sum, neg_cnt


def _weighted_bin_diff(bp_sum, bp_cnt, bn_sum, bn_cnt):
    valid = (bp_cnt > 0) & (bn_cnt > 0)
    if not valid.any():
        return np.nan
    mean_pos = bp_sum[valid] / bp_cnt[valid]
    mean_neg = bn_sum[valid] / bn_cnt[valid]
    weights = bp_cnt[valid] + bn_cnt[valid]
    return float(np.average(mean_pos - mean_neg, weights=weights))


def matched_bin_bootstrap_chunk(b_count, seed, n_prot, pos_sum, pos_cnt, neg_sum, neg_cnt):
    """Run b_count replicates of the stratified mean-diff statistic. Thread-safe
    (arrays are read-only, shared by reference)."""
    rng = np.random.default_rng(seed)
    diffs = np.empty(b_count)
    for b in range(b_count):
        c = rng.choice(n_prot, size=n_prot, replace=True)
        diffs[b] = _weighted_bin_diff(
            pos_sum[c].sum(axis=0), pos_cnt[c].sum(axis=0),
            neg_sum[c].sum(axis=0), neg_cnt[c].sum(axis=0),
        )
    return diffs


def matched_mean_diff_stats(df, prot_a_col, prot_b_col, pos_col, value_col, n_tested_col,
                             n_bins=10, B=5000, n_workers=1, seed=0):
    """
    Cluster-bootstrap CI/p-value for a log10(n_tested)-stratified mean difference
    (pos - neg), e.g. co-expression HRI vs HRNI matched on test intensity.

    Decision: bin boundaries (qcut on pooled log10(n_tested)) are computed once from
    the observed data and held fixed across all B replicates, rather than re-derived
    from each resampled multiset. This is both the standard stratified-bootstrap
    convention (strata are a property of the design, not re-estimated per replicate)
    and what makes this vectorizable: per-protein/per-bin sums are precomputed once,
    then each replicate is pure numpy array indexing - same
    precompute-once-then-resample pattern as cluster_bootstrap_stats, extended with a
    bin axis. Replaces a prior implementation that recomputed pd.qcut/groupby from
    scratch (single-threaded, O(n_edges) pandas work) on every one of B replicates.

    Returns (obs, ci_lo, ci_hi, p_val, n_bins_used).
    """
    log_n_tested = np.log10(df[n_tested_col].clip(lower=1).to_numpy())
    try:
        binned = pd.qcut(log_n_tested, q=n_bins, duplicates="drop")
    except ValueError:
        return np.nan, np.nan, np.nan, np.nan, 0
    bin_codes = binned.codes.astype(np.intp)
    n_bins_used = len(binned.categories)

    prot_a = df[prot_a_col].to_numpy()
    prot_b = df[prot_b_col].to_numpy()
    proteins = pd.unique(np.concatenate([prot_a, prot_b]))
    prot_to_i = {p: i for i, p in enumerate(proteins)}
    n_prot = len(proteins)
    prot_a_idx = np.fromiter((prot_to_i[p] for p in prot_a), dtype=np.intp, count=len(prot_a))
    prot_b_idx = np.fromiter((prot_to_i[p] for p in prot_b), dtype=np.intp, count=len(prot_b))
    pos_mask = df[pos_col].to_numpy().astype(bool)
    values = df[value_col].to_numpy(dtype=np.float64)

    pos_sum, pos_cnt, neg_sum, neg_cnt = precompute_protein_bin_sums(
        prot_a_idx, prot_b_idx, pos_mask, bin_codes, values, n_prot, n_bins_used
    )

    obs = _weighted_bin_diff(pos_sum.sum(axis=0), pos_cnt.sum(axis=0),
                              neg_sum.sum(axis=0), neg_cnt.sum(axis=0))

    q, r = divmod(B, n_workers)
    chunk_sizes = [q + (1 if i < r else 0) for i in range(n_workers)]
    chunk_seeds = [seed + i * (B + 1) for i in range(n_workers)]

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = [
            pool.submit(matched_bin_bootstrap_chunk, b_count, s, n_prot,
                        pos_sum, pos_cnt, neg_sum, neg_cnt)
            for b_count, s in zip(chunk_sizes, chunk_seeds)
        ]
        diffs = np.concatenate([f.result() for f in futures])

    diffs = diffs[~np.isnan(diffs)]
    if len(diffs) == 0 or np.isnan(obs):
        return obs, np.nan, np.nan, np.nan, n_bins_used
    ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])
    p_val = float(np.clip(2 * min((diffs < obs).mean(), (diffs > obs).mean()), 1 / len(diffs), 1))
    return obs, ci_lo, ci_hi, p_val, n_bins_used


def iid_bootstrap_ci(n_items, statistic_fn, B=5000, seed=0):
    """
    Plain (non-cluster) bootstrap for statistics computed on a table with one
    independent row per protein (no bait/prey edge dependence) - e.g. a Spearman
    correlation between per-protein interaction_ratio and n_pubmed.

    statistic_fn(idx: np.ndarray) -> float, given a resampled array of row indices
    (0..n_items-1, drawn with replacement).

    Returns (ci_lo, ci_hi, values) where values is the array of B replicate statistics.
    """
    rng = np.random.default_rng(seed)
    values = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n_items, size=n_items)
        values[b] = statistic_fn(idx)
    ci_lo, ci_hi = np.percentile(values, [2.5, 97.5])
    return ci_lo, ci_hi, values


def cluster_bootstrap_edges(df, prot_a_col, prot_b_col, statistic_fn, B=5000, seed=0):
    """
    Protein-level cluster bootstrap for edge-level statistics that are not reducible
    to per-protein sums (e.g. Cohen's kappa or an O/E ratio on a 2x2 table over pairs).
    O(n_edges) per iteration - use precompute_protein_sums/bootstrap_chunk instead
    when the statistic decomposes into per-protein sums (much faster for large B).

    statistic_fn(edge_subframe: pd.DataFrame) -> float

    Returns (obs, ci_lo, ci_hi, p_val, values).
    """
    rng = np.random.default_rng(seed)
    prot_a = df[prot_a_col].to_numpy()
    prot_b = df[prot_b_col].to_numpy()
    proteins = pd.unique(np.concatenate([prot_a, prot_b]))
    idx_by_protein = {
        p: df.index[(prot_a == p) | (prot_b == p)].to_numpy()
        for p in proteins
    }
    n_prot = len(proteins)
    obs = statistic_fn(df)

    values = np.empty(B)
    for b in range(B):
        picked = rng.choice(proteins, size=n_prot, replace=True)
        rows = np.concatenate([idx_by_protein[p] for p in picked])
        values[b] = statistic_fn(df.loc[rows])

    ci_lo, ci_hi = np.percentile(values, [2.5, 97.5])
    p_val = float(np.clip(2 * min((values < obs).mean(), (values > obs).mean()), 1 / B, 1))
    return obs, ci_lo, ci_hi, p_val, values
