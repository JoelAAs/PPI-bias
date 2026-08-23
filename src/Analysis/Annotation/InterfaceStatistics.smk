import pandas as pd



def get_interface_gene_pairs(interface_file, gene_to_id_file, id_pattern):
    interfaces = pd.read_csv(interface_file, sep="\t")
    interfaces[["uniprot_a", "uniprot_b"]] = interfaces["Protein pair"].str.split("_", expand=True)

    if id_pattern == "gene_name":
        id_map = dict()
        with open(gene_to_id_file, "r") as f:
            next(f)
            for line in f:
                uni, gn = line.strip().split("\t")
                id_map[uni] = gn
        interfaces["id_a"] = interfaces["uniprot_a"].map(id_map).fillna("NA")
        interfaces["id_b"] = interfaces["uniprot_b"].map(id_map).fillna("NA")
    else:
        interfaces["id_a"] = interfaces["uniprot_a"]
        interfaces["id_b"] = interfaces["uniprot_b"]

    interfaces["gene_pair"] = interfaces.apply(
        lambda r: "_".join(sorted([r["id_a"], r["id_b"]])), axis=1
    )
    return interfaces


def get_pair_id(row, bait_column, prey_column):
    return "_".join(sorted([row[bait_column], row[prey_column]]))


def get_mean_per_categories(pod_df, interfaces_df, bait_column, prey_column):
    """
    Per interface-size Category: the pair-wise detection ratio (mean of each
    pair's own n_observed/n_tested) vs. the global detection ratio (pooled
    sum(n_observed)/sum(n_tested)) - these can diverge when pairs within a
    category are tested very unevenly.
    """
    pod_df = pod_df.copy()
    pod_df["gene_pair"] = pod_df.apply(get_pair_id, axis=1, args=(bait_column, prey_column))
    pod_df = pod_df.merge(interfaces_df[["gene_pair", "Category"]], on="gene_pair", how="left")
    pod_df["Category"] = pod_df["Category"].fillna("Unknown")
    pod_df["pair_detection_ratio"] = pod_df["n_observed"] / pod_df["n_tested"]

    sum_stats = pod_df.groupby("Category").agg(
        pair_detection_ratio=("pair_detection_ratio", "mean"),
        total_tested_pairs=("pair_detection_ratio", "size"),
        total_observed=("n_observed", "sum"),
        total_tests=("n_tested", "sum"),
    )
    sum_stats["global_detection_ratio"] = sum_stats["total_observed"] / sum_stats["total_tests"]

    return sum_stats.reset_index()



rule get_detection_rate_per_interface_size:
    params:
        id_pattern = config["id_pattern"],
        bait_column = f"{config['id_pattern']}_bait",
        prey_column = f"{config['id_pattern']}_prey"
    input:
        interfaces = "data/DCA/benchmarks/pairs_partitioned_by_interface_sizes.tsv.gz",
        gene_to_uniprot = "work_folder/gene_names/uniprot_to_gene_name.csv",
        pod = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq"
    output:
        summary_stats = "work_folder/analysis/interfaces/{dataset}_{network_type}_detection.csv"
    log:
        "logs/analysis/interfaces/{dataset}_{network_type}_detection.log"
    run:
        interfaces_df = get_interface_gene_pairs(input.interfaces, input.gene_to_uniprot, params.id_pattern)
        pod_df = pd.read_parquet(input.pod)
        sum_stats = get_mean_per_categories(pod_df, interfaces_df, params.bait_column, params.prey_column)
        sum_stats.to_csv(output.summary_stats, sep="\t", index=False)


rule plot_interface_detection_rate:
    input:
        expand("work_folder/analysis/interfaces/{dataset}_{{network_type}}_detection.csv",
            dataset=[d for d in config["datasets"] if d != "flat"])
    output:
        "work_folder/analysis/interfaces/plots/{network_type}_interface_detection.png"
    log:
        "logs/analysis/interfaces/plots/{network_type}_interface_detection.log"
    script:
        "scripts/plot_interface_detection_rate.R"

