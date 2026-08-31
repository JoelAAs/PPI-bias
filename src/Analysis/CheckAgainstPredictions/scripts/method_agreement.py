import pandas as pd
import yaml

# load and sort
df_pod = pd.read_parquet(snakemake.input.pod,
                          columns=["uniprot_id_bait", "uniprot_id_prey", "n_tested", "n_observed"])
pod_ids = df_pod[["uniprot_id_bait", "uniprot_id_prey"]]
df_pod["pair_id"] = pod_ids.min(axis=1) + "_" + pod_ids.max(axis=1)

other_df = pd.read_csv(snakemake.input.other_ppi_intact, sep="\t")
other_ids = other_df[["prot_a", "prot_b"]]
other_df["pair_id"] = other_ids.min(axis=1) + "_" + other_ids.max(axis=1)

with open(snakemake.input.method_groupings) as f:
    method_groups = yaml.safe_load(f)

tested = df_pod["n_tested"].fillna(0)
observed = df_pod["n_observed"].fillna(-1)

interactions = df_pod[observed > snakemake.params.interaction_limit]
interaction_pairs = set(interactions["pair_id"])

hrni_sets = dict()
for k in snakemake.params.HRNI_limits:
    hrni_sets[f"HRNI_n{k}"] = df_pod[(observed == 0) & (tested >= k)]

rows = []
for group_name, mi_codes in method_groups.items():
    group_pairs = set(other_df.loc[other_df["detection_method"].isin(mi_codes), "pair_id"])
    n_group = len(group_pairs)
    interaction_overlap = len(group_pairs & interaction_pairs)

    for set_name, sub in hrni_sets.items():
        pod_pairs = set(sub["pair_id"])
        n_set = len(pod_pairs)
        n_overlap = len(pod_pairs & group_pairs)
        rows.append(dict(
            group=group_name,
            set=set_name,
            n_pod_set=n_set,
            n_other_ppi_group=n_group,
            n_overlap=n_overlap,
            discovered_disagreement=n_overlap / n_group if n_group else float("nan"),
            discovered_agreement=interaction_overlap / n_group if n_group else float("nan"),
        ))

pd.DataFrame(rows).to_csv(snakemake.output.detection_statistics, index=False)
