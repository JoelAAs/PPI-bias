rule jaccard_intra_dist_train_positive:
    input:
        edge_files = expand("work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
            permutation=range(config["n_permutations"]))
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.csv"
    log:
        "work_folder/subsets/dataset_eval/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.log"
    resources:
        mem_gb=2
    script:
        "scripts/intra_jaccard_distance.py"


rule plot_intra_permutation_jaccard_distances:
    input:
        permuted_pos = expand("work_folder/subsets/dataset_eval/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.csv",
            dataset=config["datasets"],
            neg_limit=config["negative_limits"],
            pos_limit=config["positive_limits"])    
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/jaccard_dist/{network_type}_intra_permutation_jaccard_distances.png"
    resources:
        mem_gb=2
    script:
        "scripts/plot_intra_permutation_jaccard_distance.py"


rule jaccard_inter_dist_train_positive:
    params:
        pos_limits= config["positive_limits"]
    input:
        permuted_pos = expand("work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{pos_limit}_pos.csv",
            permutation=range(config["n_permutations"]),
            pos_limit=config["positive_limits"])
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/intra_{dataset}_{network_type}_limit_{neg_limit}_poslim-jaccard_distances.csv"
    log:
        "work_folder/subsets/dataset_eval/intra_{dataset}_{network_type}_limit_{neg_limit}_poslim-jaccard_distances.log"
    resources:
        mem_gb=2
    script:
        "scripts/inter_jaccard_distance.py"


rule plot_inter_permutation_jaccard_distances:
    input:
        permuted_pos = expand("work_folder/subsets/dataset_eval/intra_{dataset}_{{network_type}}_limit_{neg_limit}_poslim-jaccard_distances.csv",
            dataset=config["datasets"],
            neg_limit=config["negative_limits"])
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/jaccard_dist/{network_type}_inter_permutation_jaccard_distances.png"
    resources:
        mem_gb=2
    script:
        "scripts/plot_inter_permutation_jaccard_distance.py"


def negative_train_edges(wildcards):
    # HRNI ("true" negatives) are the plain _neg.csv, the non-observed (NO) negatives are the
    # degree-matched -random_neg.csv drawn for the same permutation/threshold config.
    suffix = "_neg.csv" if wildcards.negative_type == "hrni" else "-random_neg.csv"
    return expand(
        "work_folder/subsets/train/permuted/{permutation}/"
        f"{wildcards.dataset}_{wildcards.network_type}_limit_{wildcards.neg_limit}"
        f"_poslim_{wildcards.pos_limit}{suffix}",
        permutation=range(config["n_permutations"]),
    )


rule jaccard_intra_dist_train_negative:
    # Same across-permutation edge overlap as jaccard_intra_dist_train_positive, but for the two
    # negative training sets (HRNI and NO) instead of the positives.
    wildcard_constraints:
        negative_type="(hrni|no)"
    input:
        edge_files = negative_train_edges
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/negatives/{negative_type}_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.csv"
    log:
        "work_folder/subsets/dataset_eval/negatives/{negative_type}_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.log"
    resources:
        mem_gb=2
    script:
        "scripts/intra_jaccard_distance.py"


rule plot_train_edge_jaccard_boxplot:
    # Three boxes, each pooling every dataset x pos_limit x neg_limit config:
    #  - HRNI negatives across the 10 permutations of one config
    #  - NO negatives across the 10 permutations of one config
    #  - positives between adjacent positive-threshold configs
    input:
        hrni_negative = expand("work_folder/subsets/dataset_eval/negatives/hrni_{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.csv",
            dataset=config["datasets"],
            neg_limit=config["negative_limits"],
            pos_limit=config["positive_limits"]),
        no_negative = expand("work_folder/subsets/dataset_eval/negatives/no_{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}-jaccard_distances.csv",
            dataset=config["datasets"],
            neg_limit=config["negative_limits"],
            pos_limit=config["positive_limits"]),
        positive_between_config = expand("work_folder/subsets/dataset_eval/intra_{dataset}_{{network_type}}_limit_{neg_limit}_poslim-jaccard_distances.csv",
            dataset=config["datasets"],
            neg_limit=config["negative_limits"])
    output:
        jaccard_distances = "work_folder/subsets/dataset_eval/jaccard_dist/{network_type}_train_edge_jaccard_boxplot.png"
    resources:
        mem_gb=2
    script:
        "scripts/plot_train_edge_jaccard_boxplot.py"
