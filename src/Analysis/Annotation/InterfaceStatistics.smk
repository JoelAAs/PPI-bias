import pandas as pd



def get_interface_gene_names_dict(interface_file, gene_name_to_uniprot):
        interfaces = pd.read_csv(interface_file, sep="\t")
        interfaces[["uniprot_a", "uniprot_b"]] = interfaces["Protein pair"].str.split("_", expand=True)
        gene_name_dict = dict()
        with open(gene_name_to_uniprot, "r") as f:
            next(f)
            for line in f:
                uni, gn = line.strip().split("\t")
                gene_name_dict[uni] = gn

        return gene_name_dict, interfaces


def get_gene_comb(row, gene_name_dict):
    gnl = [gene_name_dict.get(row["uniprot_a"], "NA"), gene_name_dict.get(row["uniprot_b"], "NA")]
    return "_".join(sorted(gnl))


def get_gene_pair(row):
    return "_".join(sorted([row["gene_name_bait"], row["gene_name_prey"]]))


def get_mean_per_categories(filename, interfaces_df):
    df_hit = pd.read_csv(filename, sep="\t")
    df_hit["gene_pair"] = df_hit.apply(get_gene_pair, axis=1)
    df_hit = df_hit.merge(interfaces_df[["gene_pair", "Category"]], on="gene_pair", how="left")
    df_hit["Category"] = df_hit["Category"].fillna("Unknown")
    df_hit["pair_detection_ratio"] = df_hit["n_observed"]/df_hit["n_tested"]
    dr_ratio = df_hit.groupby("Category")["pair_detection_ratio"].mean()
    n_pairs_tested = df_hit.groupby("Category").size()
    n_observations = df_hit.groupby("Category")["n_observed"].sum()
    n_tested = df_hit.groupby("Category")["n_tested"].sum()
    sum_stats = pd.concat([
        dr_ratio,
        n_pairs_tested,
        n_observations,
        n_tested
    ], axis=1)
    sum_stats.columns = [
        "pair_detection_ratio",
        "total_tested_pairs",
        "total_observed",
        "total_tests"
    ]

    return sum_stats.reset_index()



rule get_detection_rate_per_interface_size:
    input:
        interfaces = "data/DCA/benchmarks/pairs_partitioned_by_interface_sizes.tsv",
        gene_to_uniprot = f"work_folder{pn}/intact/uniprot_to_gene_name.csv",
        pod = f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        summary_stats = f"work_folder{pn}/analysis/interfaces/detection_{{data}}.csv"
    log:
        f"logs{pn}/analysis/interfaces/detection_{{data}}.log"
    run:
        gene_name_dict, interfaces_df = get_interface_gene_names_dict(input.interfaces, input.gene_to_uniprot)
        interfaces_df["gene_pair"] = interfaces_df[["uniprot_a", "uniprot_b"]].apply(get_gene_comb,axis=1, args=(gene_name_dict,))
        sum_stats = get_mean_per_categories(input.pod, interfaces_df)
        sum_stats.to_csv(output.summary_stats, sep="\t", index=False)

