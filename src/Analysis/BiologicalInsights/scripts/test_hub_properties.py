"""
Module 2 - hub properties (BiologicalInsights plan §7.4-7.5): derive UniProtKB-based
features for the hub set (high deg_neg, medium deg_pos) and test each for simple
over/underrepresentation against a reference population of actual interactors (low
deg_neg, medium-or-high deg_pos, same tested-depth restriction) - a plain two-sample
bootstrap, not a matched-control design.

"single_domain" (plan §7.4) requires an extra InterPro join that is not implemented
in this pass - dropped entirely rather than faked, per the plan's explicit instruction.
"""
import numpy as np
import pandas as pd

from src.Analysis.BiologicalInsights.scripts.bootstrap import benjamini_hochberg

BINARY_FEATURES = ["enzyme", "narrow_substrate_enzyme", "membrane", "secreted"]
CONTINUOUS_FEATURES = ["length", "disorder_fraction"]
EPS = 1e-10


def build_features(accessions, uniprot_df, disorder_df):
    df = uniprot_df.set_index("uniprot_id").reindex(accessions)
    df.index.name = "uniprot_id"
    df = df.reset_index()

    ec = df["ec"].fillna("")
    df["enzyme"] = ec.str.strip() != ""

    def is_narrow(row):
        ec_val = row["ec"]
        if not isinstance(ec_val, str) or not ec_val.strip():
            return False
        ec_numbers = [e.strip() for e in ec_val.split(";") if e.strip()]
        if len(ec_numbers) != 1:
            return False
        parts = ec_numbers[0].split(".")
        if len(parts) != 4 or parts[3] == "-":
            return False
        cat = row.get("cc_catalytic_activity")
        n_activities = str(cat).count("CATALYTIC ACTIVITY:") if isinstance(cat, str) else 0
        return n_activities == 1

    df["narrow_substrate_enzyme"] = df.apply(is_narrow, axis=1)
    df["membrane"] = df["ft_transmem"].fillna("").str.strip() != ""
    subloc = df["cc_subcellular_location"].fillna("")
    signal = df["ft_signal"].fillna("")
    df["secreted"] = (signal.str.strip() != "") | subloc.str.contains("Secreted", na=False)

    disorder_map = disorder_df.set_index("uniprot_id")["disorder_fraction"]
    df["disorder_fraction"] = df["uniprot_id"].map(disorder_map)
    return df


def two_sample_bootstrap_or(hub_vals, pool_vals, B, seed):
    rng = np.random.default_rng(seed)
    r_hub = hub_vals.mean()
    r_pool = pool_vals.mean()
    obs_or = (r_hub / max(1 - r_hub, EPS)) / max(r_pool / max(1 - r_pool, EPS), EPS)
    n_hub, n_pool = len(hub_vals), len(pool_vals)
    log_ors = np.empty(B)
    for b in range(B):
        rh = hub_vals[rng.integers(0, n_hub, n_hub)].mean()
        rp = pool_vals[rng.integers(0, n_pool, n_pool)].mean()
        oh = rh / max(1 - rh, EPS)
        op = max(rp / max(1 - rp, EPS), EPS)
        log_ors[b] = np.log(max(oh, EPS) / op)
    ci_lo, ci_hi = np.exp(np.percentile(log_ors, [2.5, 97.5]))
    p = float(np.clip(2 * min((log_ors < 0).mean(), (log_ors > 0).mean()), 1 / B, 1))
    return float(obs_or), float(ci_lo), float(ci_hi), p


def two_sample_bootstrap_pct_diff(hub_vals, pool_vals, B, seed):
    """% difference in the mean (hub vs reference), relative to the reference mean."""
    rng = np.random.default_rng(seed)
    r_hub = hub_vals.mean()
    r_pool = pool_vals.mean()
    obs_pct = float((r_hub - r_pool) / max(abs(r_pool), EPS) * 100)
    n_hub, n_pool = len(hub_vals), len(pool_vals)
    pct_diffs = np.empty(B)
    for b in range(B):
        rh = hub_vals[rng.integers(0, n_hub, n_hub)].mean()
        rp = pool_vals[rng.integers(0, n_pool, n_pool)].mean()
        pct_diffs[b] = (rh - rp) / max(abs(rp), EPS) * 100
    ci_lo, ci_hi = np.percentile(pct_diffs, [2.5, 97.5])
    p = float(np.clip(2 * min((pct_diffs < 0).mean(), (pct_diffs > 0).mean()), 1 / B, 1))
    return obs_pct, float(ci_lo), float(ci_hi), p


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    B = snakemake.params.bootstrap_replicates
    seed = snakemake.params.seed

    hub_set = pd.read_csv(snakemake.input.hub_set, sep="\t")
    reference_set = pd.read_csv(snakemake.input.reference_set, sep="\t")
    uniprot_df = pd.read_csv(snakemake.input.uniprot_annotation, sep="\t")
    disorder_df = pd.read_csv(snakemake.input.disorder, sep="\t")

    if len(hub_set) == 0 or len(reference_set) == 0:
        print("Hub set or reference set is empty - writing an empty "
              "enrichment table, no testing performed.", file=log, flush=True)
        pd.DataFrame(columns=["feature", "feature_type", "effect", "ci_lo", "ci_hi",
                               "p_val", "q_val", "n_hub", "n_reference"]).to_csv(
            snakemake.output.enrichment, sep="\t", index=False
        )
        log.close()
        raise SystemExit(0)

    hub_ids = hub_set["uniprot_id"].tolist()
    reference_ids = reference_set["uniprot_id"].tolist()
    all_ids = sorted(set(hub_ids) | set(reference_ids))
    print(f"{len(hub_ids)} hubs, {len(reference_ids)} reference (actual-interactor) proteins",
          file=log, flush=True)

    features_df = build_features(all_ids, uniprot_df, disorder_df).set_index("uniprot_id")
    print("## Decisions\n"
          "'single_domain' (plan §7.4) requires an extra InterPro join not "
          "implemented here - dropped entirely, not faked with a placeholder.\n"
          "Reference population is deg_neg_bin == 'low' AND deg_pos_bin in ('medium', "
          "'high') - actual interactors, same tested-depth restriction as the hub set "
          "- tested with a plain two-sample bootstrap, not matched controls.",
          file=log, flush=True)

    rows = []
    for feat in BINARY_FEATURES:
        hub_vals = features_df.loc[hub_ids, feat].to_numpy(dtype=float)
        ref_vals = features_df.loc[reference_ids, feat].to_numpy(dtype=float)
        obs_or, ci_lo, ci_hi, p_val = two_sample_bootstrap_or(hub_vals, ref_vals, B, seed)
        print(f"{feat}: OR={obs_or:.3f} [{ci_lo:.3f},{ci_hi:.3f}] p={p_val:.4g}",
              file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "binary", "effect": obs_or,
                      "ci_lo": ci_lo, "ci_hi": ci_hi, "p_val": p_val,
                      "n_hub": len(hub_ids), "n_reference": len(reference_ids)})

    for feat in CONTINUOUS_FEATURES:
        hub_vals = features_df.loc[hub_ids, feat].dropna().to_numpy(dtype=float)
        ref_vals = features_df.loc[reference_ids, feat].dropna().to_numpy(dtype=float)
        obs_pct, ci_lo, ci_hi, p_val = two_sample_bootstrap_pct_diff(hub_vals, ref_vals, B, seed)
        print(f"{feat}: hub vs reference % diff={obs_pct:.4g}% [{ci_lo:.4g},{ci_hi:.4g}] "
              f"p={p_val:.4g}", file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "continuous", "effect": obs_pct,
                      "ci_lo": ci_lo, "ci_hi": ci_hi, "p_val": p_val,
                      "n_hub": len(hub_vals), "n_reference": len(ref_vals)})

    result_df = pd.DataFrame(rows)
    result_df["q_val"] = benjamini_hochberg(result_df["p_val"].to_numpy())

    print("## Interpretation (stated before looking at effect directions)\n"
          "Prior: enriched for specific enzymes, smaller/single-domain, membrane and "
          "secreted, depleted in disorder - OR membrane could instead be depleted "
          "because of detection difficulty. Either direction is a reportable result.",
          file=log, flush=True)

    result_df.to_csv(snakemake.output.enrichment, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
