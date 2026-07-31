"""
Module 2 - hub properties (BiologicalInsights plan §7.4-7.5): derive UniProtKB-based
features and test each for over/underrepresentation in non-interaction hubs
(hub_status == "non-interaction") vs actual interactors (hub_status == "interaction").
Fisher's exact test for binary features, Mann-Whitney U for continuous - both give an
exact/asymptotic p-value directly, no resampling needed.

"single_domain" (plan §7.4) requires an extra InterPro join that is not implemented
in this pass - dropped entirely rather than faked, per the plan's explicit instruction.
"""
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from src.Analysis.BiologicalInsights.scripts.bootstrap import benjamini_hochberg

BINARY_FEATURES = ["enzyme", "narrow_substrate_enzyme", "membrane", "secreted"]
CONTINUOUS_FEATURES = ["length", "disorder_fraction"]


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


def fisher_or(hub_vals, ref_vals):
    table = [[hub_vals.sum(), (~hub_vals).sum()], [ref_vals.sum(), (~ref_vals).sum()]]
    obs_or, p_val = fisher_exact(table)
    return float(obs_or), float(p_val)


def mannwhitney_pct_diff(hub_vals, ref_vals):
    """% difference in the mean (hub vs reference), relative to the reference mean."""
    obs_pct = float((hub_vals.mean() - ref_vals.mean()) / abs(ref_vals.mean()) * 100)
    _, p_val = mannwhitneyu(hub_vals, ref_vals, alternative="two-sided")
    return obs_pct, float(p_val)


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    hub_set = pd.read_csv(snakemake.input.hub_set, sep="\t")
    uniprot_df = pd.read_csv(snakemake.input.uniprot_annotation, sep="\t")
    disorder_df = pd.read_csv(snakemake.input.disorder, sep="\t")

    hub_ids = hub_set.loc[hub_set["hub_status"] == "non-interaction", "uniprot_id"].tolist()
    reference_ids = hub_set.loc[hub_set["hub_status"] == "interaction", "uniprot_id"].tolist()
    print(f"{len(hub_ids)} non-interaction hubs, {len(reference_ids)} interaction "
          f"reference proteins", file=log, flush=True)

    if not hub_ids or not reference_ids:
        print("Hub set or reference set is empty - writing an empty "
              "enrichment table, no testing performed.", file=log, flush=True)
        pd.DataFrame(columns=["feature", "feature_type", "effect", "p_val", "q_val",
                               "n_hub", "n_reference"]).to_csv(
            snakemake.output.enrichment, sep="\t", index=False
        )
        log.close()
        raise SystemExit(0)

    all_ids = sorted(set(hub_ids) | set(reference_ids))
    features_df = build_features(all_ids, uniprot_df, disorder_df).set_index("uniprot_id")

    rows = []
    for feat in BINARY_FEATURES:
        hub_vals = features_df.loc[hub_ids, feat].to_numpy(dtype=bool)
        ref_vals = features_df.loc[reference_ids, feat].to_numpy(dtype=bool)
        obs_or, p_val = fisher_or(hub_vals, ref_vals)
        print(f"{feat}: OR={obs_or:.3f} p={p_val:.4g}", file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "binary", "effect": obs_or,
                      "p_val": p_val, "n_hub": len(hub_ids), "n_reference": len(reference_ids)})

    for feat in CONTINUOUS_FEATURES:
        hub_vals = features_df.loc[hub_ids, feat].dropna().to_numpy(dtype=float)
        ref_vals = features_df.loc[reference_ids, feat].dropna().to_numpy(dtype=float)
        obs_pct, p_val = mannwhitney_pct_diff(hub_vals, ref_vals)
        print(f"{feat}: hub vs reference % diff={obs_pct:.4g}% p={p_val:.4g}",
              file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "continuous", "effect": obs_pct,
                      "p_val": p_val, "n_hub": len(hub_vals), "n_reference": len(ref_vals)})

    result_df = pd.DataFrame(rows)
    result_df["q_val"] = benjamini_hochberg(result_df["p_val"].to_numpy())

    result_df.to_csv(snakemake.output.enrichment, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
