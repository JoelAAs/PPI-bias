import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu


BINARY_FEATURES = ["enzyme", "narrow_substrate_enzyme", "membrane", "secreted"]
CONTINUOUS_FEATURES = ["length", "disorder_fraction"]

def benjamini_hochberg(p_vals):
    p_vals = np.asarray(p_vals)
    n = len(p_vals)
    order = np.argsort(p_vals)
    ranks_order = np.arange(1, n + 1)
    q_vals = np.minimum.accumulate((p_vals[order] * n / ranks_order)[::-1])[::-1][np.argsort(order)]
    return np.clip(q_vals, 0, 1)



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


def fisher_or(ni_hub_vals, i_hub_vals):
    table = [[ni_hub_vals.sum(), (~ni_hub_vals).sum()], [i_hub_vals.sum(), (~i_hub_vals).sum()]]
    obs_or, p_val = fisher_exact(table)
    return float(obs_or), float(p_val)


def mannwhitney_pct_diff(ni_hub_vals, i_hub_vals):
    """% difference in the mean (ni_hub vs i_hub), relative to the i_hub mean."""
    obs_pct = float((ni_hub_vals.mean() - i_hub_vals.mean()) / abs(i_hub_vals.mean()) * 100)
    _, p_val = mannwhitneyu(ni_hub_vals, i_hub_vals, alternative="two-sided")
    return obs_pct, float(p_val)


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    hub_set = pd.read_csv(snakemake.input.hub_set, sep="\t")
    uniprot_df = pd.read_csv(snakemake.input.uniprot_annotation, sep="\t")
    disorder_df = pd.read_csv(snakemake.input.disorder, sep="\t")

    ni_hub_ids = hub_set.loc[hub_set["hub_status"] == "non-interaction", "uniprot_id"].tolist()
    i_hub_ids = hub_set.loc[hub_set["hub_status"] == "interaction", "uniprot_id"].tolist()
    print(f"{len(ni_hub_ids)} non-interaction hubs (ni_hub), {len(i_hub_ids)} interaction "
          f"hubs (i_hub)", file=log, flush=True)

    if not ni_hub_ids or not i_hub_ids:
        print("ni_hub set or i_hub set is empty - writing an empty "
              "enrichment table, no testing performed.", file=log, flush=True)
        pd.DataFrame(columns=["feature", "feature_type", "effect", "p_val", "q_val",
                               "n_ni_hub", "n_i_hub"]).to_csv(
            snakemake.output.enrichment, sep="\t", index=False
        )
        log.close()
        raise SystemExit(0)

    all_ids = sorted(set(ni_hub_ids) | set(i_hub_ids))
    features_df = build_features(all_ids, uniprot_df, disorder_df).set_index("uniprot_id")

    rows = []
    for feat in BINARY_FEATURES:
        ni_hub_vals = features_df.loc[ni_hub_ids, feat].to_numpy(dtype=bool)
        i_hub_vals = features_df.loc[i_hub_ids, feat].to_numpy(dtype=bool)
        obs_or, p_val = fisher_or(ni_hub_vals, i_hub_vals)
        print(f"{feat}: OR={obs_or:.3f} p={p_val:.4g}", file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "binary", "effect": obs_or,
                      "p_val": p_val, "n_ni_hub": len(ni_hub_ids), "n_i_hub": len(i_hub_ids)})

    for feat in CONTINUOUS_FEATURES:
        ni_hub_vals = features_df.loc[ni_hub_ids, feat].dropna().to_numpy(dtype=float)
        i_hub_vals = features_df.loc[i_hub_ids, feat].dropna().to_numpy(dtype=float)
        obs_pct, p_val = mannwhitney_pct_diff(ni_hub_vals, i_hub_vals)
        print(f"{feat}: ni_hub vs i_hub % diff={obs_pct:.4g}% p={p_val:.4g}",
              file=log, flush=True)
        rows.append({"feature": feat, "feature_type": "continuous", "effect": obs_pct,
                      "p_val": p_val, "n_ni_hub": len(ni_hub_vals), "n_i_hub": len(i_hub_vals)})

    result_df = pd.DataFrame(rows)
    result_df["q_val"] = benjamini_hochberg(result_df["p_val"].to_numpy())

    result_df.to_csv(snakemake.output.enrichment, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
