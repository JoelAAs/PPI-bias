import pandas as pd

gene_name_df = pd.read_csv(f"work_folder/uniprot_pod/intact/uniprot_to_gene_name.csv", sep="\t")


df_ms = pd.read_csv("work_folder/uniprot_pod/analysis/POD/POD_ms.csv", sep="\t")
df_ms = df_ms.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
del df_ms["uniprot_id"]
df_ms = df_ms.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id",
                                  suffixes=("_bait", "_prey"))
del df_ms["uniprot_id"]

n_ids_baits = df_ms.groupby("gene_name_bait", as_index=False).nunique("uniprot_id_bait")
n_ids_prey = df_ms.groupby("gene_name_prey", as_index=False).nunique("uniprot_id_prey")

print(f'The number of baits with more than 1 uniprot identifiers per gene:')
print(f'{sum(n_ids_baits["size"] != 1)} / {sum(n_ids_baits["size"] != 1)}')

print(f'The number of prey with more than 1 uniprot identifiers per gene:')
print(f'{sum(n_ids_prey["size"] != 1)} / {sum(n_ids_prey["size"] != 1)}')

bait_multi_ids = n_ids_baits[n_ids_baits["size"] != 1, "gene_name_bait"]
df_ms_multi_id = df_ms[df_ms["gene_name_bait"].isin(bait_multi_ids)]

def number_of_isoforms(id_list):
    isoforms = {}
    for id in id_list:
        if "-" in id:
            isoforms &= {id.split("-")[1],}

    return len(isoforms), len(id_list) - len(isoforms),

bait_statistics = pd.DataFrame(
    gene_name = bait_multi_ids
)
bait_statistics[["n_isoforms", "n_identifiers"]] = bait_multi_ids["gene_name_bait"].apply(
    lambda x: number_of_isoforms(df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["uniprot_id"]))
bait_statistics["n_tests"] =  bait_multi_ids["gene_name_bait"].apply(
    lambda x: (df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["n_tested"].sum()))
bait_statistics["n_observed"] =  bait_multi_ids["gene_name_bait"].apply(
    lambda x: (df_ms_multi_id[df_ms_multi_id["gene_name_bait"] == x]["n_observed"].sum()))
