import pandas as pd

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

bait_multi_ids = n_ids_baits[n_ids_baits["uniprot_id_bait"] != 1]["gene_name_bait"]
df_ms_multi_id = df_ms[df_ms["gene_name_bait"].isin(bait_multi_ids)]

df_ms_multi_id_dedup = df_ms_multi_id[~df_ms_multi_id["uniprot_id_bait"].duplicated()]
df_ms_multi_id_dedup["isoform"] = df_ms_multi_id_dedup["uniprot_id_bait"].apply(lambda x: "-" in x)
n_isoforms_bait = df_ms_multi_id_dedup.groupby("gene_name_bait")["isoform"].sum()
all_entries = df_ms_multi_id_dedup.groupby("gene_name_bait").size()
n_interactions = df_ms_multi_id.groupby("gene_name_bait").size().rename({"size": "n_pairs_tested"})

bait_statistics = n_isoforms_bait.join(all_entries)
bait_statistics.columns = ["n_isoform_labels", "n_uniprot_id"]
bait_statistics = bait_statistics.join(n_interactions)

bait_statistics = pd.DataFrame(
    bait_multi_ids
)
bait_statistics[["n_isoforms", "n_identifiers"]] = bait_multi_ids.apply(
    lambda x: number_of_isoforms(df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["uniprot_id_bait"]))
bait_statistics["n_tests"] = bait_multi_ids.apply(
    lambda x: (df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["n_tested"].sum()))
bait_statistics["n_observed"] = bait_multi_ids.apply(
    lambda x: (df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["n_observed"].sum()))
