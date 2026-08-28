rule jaccard_intra_dist_train_positive:
    input:
        permuted_pos = expand("work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
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
