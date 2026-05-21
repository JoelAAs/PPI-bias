import pandas as pd
import itertools
import random
import networkx as nx


# rule get_directional_balance_report:
#     input:
#         set_pos=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
#         set_neg=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
#         train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
#         validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
#         test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
#         train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
#         validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
#         test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
#     output:
#         edge_statistics = f"work_folder{pn}/subsets/balance_data/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_stats.csv"
#     log:
#         f"logs{pn}/subsets/balance_data/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_stats.log"
#     script:
#         "scripts/get_degree_metrics.py"

rule generate_test_validation_test:
    input:
        interaction_data = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        max_negative = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv"
    output:
        validation_pos = f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_pos.csv",
        validation_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_neg.csv",
        test_pos = f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_pos.csv",
        test_neg = f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_neg.csv",
    log:
        f"logs{pn}/subsets/validation/{{dataset}}_{{network_type}}.log"
    script:
        "scripts/define_validation_test.py"


rule define_max_sets:
    params:
        max_positive = config["positive_max"],
        max_negative = config["negative_max"],
    input:
        detection_df = f"work_folder{pn}/analysis/POD/{{network_type}}/POD_{{dataset}}.pq"
    output:
        max_positive = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        max_negative = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv"
    log:
        f"logs{pn}/subsets/{{dataset}}_{{network_type}}_max_sets.log"
    run:
        detection_df = pd.read_parquet(input.detection_df)
        
        positive_edges = detection_df[detection_df["lower_bound_pod"] >= params.max_positive]
        negative_edges = detection_df[(detection_df["n_tested"] >= params.max_negative) & (detection_df["n_observed"] == 0)]
        positive_edges.to_csv(output.max_positive, index=False, columns=["gene_name_bait", "gene_name_prey"], sep="\t")
        negative_edges.to_csv(output.max_negative, index=False, columns=["gene_name_bait", "gene_name_prey"], sep="\t")


rule balance_directional_to_equal_samples:
    params:
        positive_limits = config["positive_limits"],
        negative_limits = config["negative_limits"]
    input:
        test_set = f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_pos.csv",
        validation_set = f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_pos.csv",
        full_detection=f"work_folder{pn}/analysis/POD/{{network_type}}/POD_{{dataset}}.pq"
    threads: 10
    resources:
        mem_gb=50
    output:
        balanced_edges_positive = expand(
            "work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
            pn=pn, pos_limit=config["positive_limits"], neg_limit=config["negative_limits"]),
        balanced_edges_negative = expand(
            "work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}_neg.csv",
            pn=pn, pos_limit=config["positive_limits"], neg_limit=config["negative_limits"])
    log:
        f"logs{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_balance.log"
    script:
        "scripts/sample_balance_multi_network.py"



rule generate_directional_permuted_random_negative:
    input:
        balanced_positive_edges = f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
    output:
        random_negative_edges = f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}-random_neg.csv"
    log:
        f"logs{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_random_neg.log"
    resources:
        mem_gb=30
    script:
        "scripts/random_negative.py" # TODO check


rule permute_directional_balanced_samples:
    params:
        fraction  = config.get("permutation_fraction", 0.95),
        base_seed = config["seed"],
        min_flow  = config.get("permutation_min_flow", 0.95)
    input:
        balanced_pos = f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        balanced_neg = f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    output:
        permuted_pos = f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        permuted_neg = f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    resources:
        mem_gb=30
    log:
        f"logs{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.log"
    script:
        "scripts/permute_balanced_set.py" #TODO CHECK


rule generate_directional_negative_sample:
    input:
        balanced_positive_edges = f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
    output:
        random_negative_edges = f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}-random_neg.csv"
    log:
        f"logs{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_random_neg.log"
    resources:
        mem_gb=45
    script:
        "scripts/random_negative.py"