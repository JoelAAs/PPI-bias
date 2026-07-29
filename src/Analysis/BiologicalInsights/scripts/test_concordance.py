"""
Module 4 - is agreement between MS and Y2H concentrated on negatives?
Rationale (carried into the write-up): weak/transient and co-complex interactions
are both real, so positives are where the two assays differ in sensitivity, whereas
a genuine non-interaction has nothing for either assay to detect - agreement should
therefore concentrate on negatives.

{positive_definition} wildcard:
  any_pos - n_observed >= 1 in that assay
  strict  - lower_bound_pod > config["positive_max"] in that assay

Protein-level bootstrap here uses a different resampling rule than
bootstrap.cluster_bootstrap_edges, per the plan's explicit instruction for this
module: draw proteins with replacement, then keep a pair only if BOTH endpoints are
in the drawn set (no edge-multiplicity weighting). Implemented locally as
vertex_subsample_bootstrap since it is not reused elsewhere.
"""
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

N_BOOTSTRAP_MAX_EDGES_WARN = 2_000_000  # vertex_subsample_bootstrap is O(n_edges) per iter


def _bootstrap_chunk(df, prot_a_col, prot_b_col, statistic_fn, proteins, b_count, seed):
    """
    One worker's share of the B replicates. The per-iteration mask construction below
    is a pure-Python loop (holds the GIL), so this is parallelised with
    ProcessPoolExecutor, not threads - unlike bootstrap.py's numpy-vectorised chunks,
    threads here would not actually run concurrently.
    """
    rng = np.random.default_rng(seed)
    n_prot = len(proteins)
    prot_a = df[prot_a_col].to_numpy()
    prot_b = df[prot_b_col].to_numpy()
    values = np.empty(b_count)
    for b in range(b_count):
        drawn = set(rng.choice(proteins, size=n_prot, replace=True))
        mask = np.fromiter((a in drawn and c in drawn for a, c in zip(prot_a, prot_b)),
                            dtype=bool, count=len(df))
        values[b] = statistic_fn(df.loc[mask])
    return values


def vertex_subsample_bootstrap(df, prot_a_col, prot_b_col, statistic_fn, B, seed, n_workers=1):
    proteins = pd.unique(np.concatenate([df[prot_a_col].to_numpy(), df[prot_b_col].to_numpy()]))
    obs = statistic_fn(df)

    if n_workers <= 1:
        values = _bootstrap_chunk(df, prot_a_col, prot_b_col, statistic_fn, proteins, B, seed)
    else:
        q, r = divmod(B, n_workers)
        chunk_sizes = [q + (1 if i < r else 0) for i in range(n_workers)]
        chunk_seeds = [seed + i * (B + 1) for i in range(n_workers)]
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = [
                pool.submit(_bootstrap_chunk, df, prot_a_col, prot_b_col, statistic_fn,
                            proteins, b_count, s)
                for b_count, s in zip(chunk_sizes, chunk_seeds)
            ]
            values = np.concatenate([f.result() for f in futures])

    ci_lo, ci_hi = np.percentile(values, [2.5, 97.5])
    return obs, ci_lo, ci_hi, values


def oe_neg(df):
    p_neg_ms = (~df["ms_pos"]).mean()
    p_neg_y2h = (~df["y2h_pos"]).mean()
    expected = p_neg_ms * p_neg_y2h * len(df)
    observed = (~df["ms_pos"] & ~df["y2h_pos"]).sum()
    return observed / expected if expected > 0 else np.nan


def oe_pos(df):
    p_pos_ms = df["ms_pos"].mean()
    p_pos_y2h = df["y2h_pos"].mean()
    expected = p_pos_ms * p_pos_y2h * len(df)
    observed = (df["ms_pos"] & df["y2h_pos"]).sum()
    return observed / expected if expected > 0 else np.nan


def cohens_kappa(df):
    n = len(df)
    p_pos_ms = df["ms_pos"].mean()
    p_pos_y2h = df["y2h_pos"].mean()
    po = ((df["ms_pos"] == df["y2h_pos"])).mean()
    pe = p_pos_ms * p_pos_y2h + (1 - p_pos_ms) * (1 - p_pos_y2h)
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


def match_prevalence(df, rng):
    """Randomly remove positive-calling pairs from the higher-prevalence assay until
    both assays' marginal positive rate equal p* = min(p_ms, p_y2h)."""
    p_ms = df["ms_pos"].mean()
    p_y2h = df["y2h_pos"].mean()
    if np.isclose(p_ms, p_y2h):
        return df
    higher_col, target_p = ("ms_pos", p_y2h) if p_ms > p_y2h else ("y2h_pos", p_ms)
    is_pos = df[higher_col].to_numpy()
    n_pos, n_neg = is_pos.sum(), (~is_pos).sum()
    if n_pos == 0:
        return df
    new_n_pos = min(n_pos, int(round(target_p * n_neg / (1 - target_p)))) # Expected number of postives if MS were as trigger happy as Y2H
    pos_idx = np.flatnonzero(is_pos)
    keep_pos_idx = rng.choice(pos_idx, size=new_n_pos, replace=False) # get a similar representative sample
    neg_idx = np.flatnonzero(~is_pos)
    keep_idx = np.concatenate([keep_pos_idx, neg_idx])
    return df.iloc[keep_idx]


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    config = snakemake.params.config
    positive_definition = snakemake.wildcards.positive_definition
    seed = config["seed"]
    B = config["bootstrap_replicates"]
    n_draws = config["n_matched_draws"]
    rule_seed = (seed + (hash(positive_definition) % 10_000)) % (2**31 - 1)
    rng = np.random.default_rng(rule_seed)

    joint = pd.read_parquet(snakemake.input.joint_space)
    print(f"{len(joint)} pairs in the joint (both-assay-tested) space, "
          f"positive_definition={positive_definition}", file=log, flush=True)

    if positive_definition == "any_pos":
        joint["ms_pos"] = joint["n_observed_ms"] >= 1
        joint["y2h_pos"] = joint["n_observed_y2h"] >= 1
    elif positive_definition == "strict":
        joint["ms_pos"] = joint["lower_bound_pod_ms"] > config["positive_max"]
        joint["y2h_pos"] = joint["lower_bound_pod_y2h"] > config["positive_max"]
    else:
        raise ValueError(f"unknown positive_definition wildcard: {positive_definition!r}")

    n_both_pos = int((joint["ms_pos"] & joint["y2h_pos"]).sum())
    n_both_neg = int((~joint["ms_pos"] & ~joint["y2h_pos"]).sum())
    n_ms_only = int((joint["ms_pos"] & ~joint["y2h_pos"]).sum())
    n_y2h_only = int((~joint["ms_pos"] & joint["y2h_pos"]).sum())
    print(f"2x2: both_pos={n_both_pos}, both_neg={n_both_neg}, "
          f"ms_only_pos={n_ms_only}, y2h_only_pos={n_y2h_only}", file=log, flush=True)
    if len(joint) > N_BOOTSTRAP_MAX_EDGES_WARN:
        print(f"WARNING: {len(joint)} edges exceeds {N_BOOTSTRAP_MAX_EDGES_WARN} - "
              f"vertex_subsample_bootstrap is O(n_edges) per iteration and B={B} "
              f"replicates may be slow", file=log, flush=True)

    rows = []
    for name, stat_fn in [("oe_negative", oe_neg), ("oe_positive", oe_pos), ("kappa", cohens_kappa)]:
        obs, ci_lo, ci_hi, _ = vertex_subsample_bootstrap(
            joint, "prot_a", "prot_b", stat_fn, B=B, seed=rule_seed, n_workers=snakemake.threads
        )
        rows.append({"statistic": name, "estimate_type": "raw", "value": obs,
                      "ci_lo": ci_lo, "ci_hi": ci_hi, "n_pairs": len(joint)})
    print("Raw (unmatched) O/E and kappa computed", file=log, flush=True)

    # prevalence fix: subsample to p* = min(p_ms, p_y2h), n_matched_draws times ----
    p_ms, p_y2h = joint["ms_pos"].mean(), joint["y2h_pos"].mean()
    print(f"Raw positive rates: p_ms={p_ms:.4g}, p_y2h={p_y2h:.4g} - matching to "
          f"p*={min(p_ms, p_y2h):.4g}", file=log, flush=True)

    draw_values = {"oe_negative": [], "oe_positive": [], "kappa": []}
    for d in range(n_draws):
        matched = match_prevalence(joint, rng)
        draw_values["oe_negative"].append(oe_neg(matched))
        draw_values["oe_positive"].append(oe_pos(matched))
        draw_values["kappa"].append(cohens_kappa(matched))

    for name in ["oe_negative", "oe_positive", "kappa"]:
        vals = np.asarray(draw_values[name])
        ci_lo, ci_hi = np.percentile(vals, [2.5, 97.5])
        rows.append({"statistic": name, "estimate_type": "prevalence_matched",
                      "value": float(np.mean(vals)), "ci_lo": ci_lo, "ci_hi": ci_hi,
                      "n_pairs": len(joint)})
    print(f"Prevalence-matched estimates computed over {n_draws} draws "
          f"(headline number per plan §9.3)", file=log, flush=True)

    result_df = pd.DataFrame(rows)
    result_df["n_both_pos"] = n_both_pos
    result_df["n_both_neg"] = n_both_neg
    result_df["n_ms_only_pos"] = n_ms_only
    result_df["n_y2h_only_pos"] = n_y2h_only
    result_df["p_ms"] = p_ms
    result_df["p_y2h"] = p_y2h
    result_df.to_csv(snakemake.output.concordance, sep="\t", index=False)

    print("Done.", file=log, flush=True)
    log.close()
