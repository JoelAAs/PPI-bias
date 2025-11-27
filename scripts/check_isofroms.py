import pandas as pd
"""
Script checking the number of uniprot identifiers and isoforms used within the MS ppi studies
"""
gene_name_df = pd.read_csv(f"work_folder/uniprot_pod/intact/uniprot_to_gene_name.csv", sep="\t")

df_ms = pd.read_csv("work_folder/uniprot_pod/analysis/POD/POD_ms.csv", sep="\t")
df_ms = df_ms.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
del df_ms["uniprot_id"]
df_ms = df_ms.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id",
                    suffixes=("_bait", "_prey"))
del df_ms["uniprot_id"]

n_ids_baits = df_ms.groupby("gene_name_bait", as_index=False)["uniprot_id_bait"].nunique()
n_ids_prey = df_ms.groupby("gene_name_prey", as_index=False)["uniprot_id_prey"].nunique()

print(f'The number of baits with more than 1 uniprot identifiers per gene:')
print(f'{sum(n_ids_baits["uniprot_id_bait"] != 1)} / {len(n_ids_baits["uniprot_id_bait"])}')

print(f'The number of prey with more than 1 uniprot identifiers per gene:')
print(f'{sum(n_ids_prey["uniprot_id_prey"] != 1)} / {len(n_ids_prey["uniprot_id_prey"])}')

# Bait
bait_multi_ids = n_ids_baits[n_ids_baits["uniprot_id_bait"] != 1]["gene_name_bait"]
df_ms_multi_id = df_ms[df_ms["gene_name_bait"].isin(bait_multi_ids)]

df_ms_multi_id_dedup = df_ms_multi_id[~df_ms_multi_id["uniprot_id_bait"].duplicated()]
df_ms_multi_id_dedup["isoform"] = df_ms_multi_id_dedup["uniprot_id_bait"].apply(lambda x: "-" in x)
n_isoforms_bait = df_ms_multi_id_dedup.groupby("gene_name_bait")["isoform"].sum()
all_entries = df_ms_multi_id_dedup.groupby("gene_name_bait").size()
n_interactions = df_ms_multi_id.groupby("gene_name_bait").size()
n_observations = df_ms_multi_id.groupby("gene_name_bait")["n_observed"].sum()
n_tests = df_ms_multi_id.groupby("gene_name_bait")["n_tested"].sum()


bait_statistics = pd.concat([
    n_isoforms_bait,
    all_entries,
    n_interactions,
    n_observations,
    n_tests
], axis = 1)
bait_statistics.columns = ["n_isoform_labels", "n_uniprot_id", "n_pairs_tested", "n_observations", "n_tests"]
bait_statistics.reset_index().to_csv("work_folder/uniprot_pod/isoforms_bait.csv")

# Prey
prey_multi_ids = n_ids_prey[n_ids_prey["uniprot_id_prey"] != 1]["gene_name_prey"]
df_ms_multi_id = df_ms[df_ms["gene_name_prey"].isin(prey_multi_ids)]

df_ms_multi_id_dedup = df_ms_multi_id[~df_ms_multi_id["uniprot_id_prey"].duplicated()]
df_ms_multi_id_dedup["isoform"] = df_ms_multi_id_dedup["uniprot_id_prey"].apply(lambda x: "-" in x)
n_isoforms_prey = df_ms_multi_id_dedup.groupby("gene_name_prey")["isoform"].sum()
all_entries = df_ms_multi_id_dedup.groupby("gene_name_prey").size()
n_interactions = df_ms_multi_id.groupby("gene_name_prey").size()
n_observations = df_ms_multi_id.groupby("gene_name_prey")["n_observed"].sum()
n_tests = df_ms_multi_id.groupby("gene_name_prey")["n_tested"].sum()


prey_statistics = pd.concat([
    n_isoforms_prey,
    all_entries,
    n_interactions,
    n_observations,
    n_tests
], axis = 1)
prey_statistics.columns = ["n_isoform_labels", "n_uniprot_id", "n_pairs_tested", "n_observations", "n_tests"]
prey_statistics.reset_index().to_csv("work_folder/uniprot_pod/isoforms_prey.csv")
