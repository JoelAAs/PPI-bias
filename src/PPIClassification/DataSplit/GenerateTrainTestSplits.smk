import pandas as pd
import itertools
import random
import networkx as nx


def get_gene_partition(wc):
    if wc["partition_name"] == "sequencesimilarity":
        return f"work_folder{pn}/subsets/partitions/sequencesimilarity_gene_name.txt"
    elif wc["partition_name"] == "maxpos":
        return f"work_folder{pn}/subsets/partitions/{wc['dataset']}_{wc["network_type"]}_limit_{wc['pos_limit']}_gene_name.txt"
    elif wc["partition_name"] == "pre-balanced":
        return f"work_folder{pn}/subsets/partitions/{wc['dataset']}__directional_limit_{wc['neg_limit']}_poslim_{wc['pos_limit']}_gene_name.txt"
    else:
        raise ValueError(f"unknown partition name = {wc['partition_name']}")


rule get_directional_splits:
    input:
        set_pos=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
    output:
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    script:
        "scripts/split_greedy.py"


rule get_directional_balance_report:
    input:
        set_pos=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    output:
        edge_statistics = f"work_folder{pn}/subsets/balance_data/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_stats.csv"
    script:
        "scripts/get_degree_metrics.py"

rule generate_validation_test:
    input:
        interaction_data = f"work_folder{pn}/subsets/{{dataset}}_directional_limit_{config['positive_max']}_pos.csv",
        max_negative = f"work_folder{pn}/subsets/{{dataset}}_directional_limit_{config['negative_max']}_pos.csv"
    output:
        validation_pos = f"work_folder{pn}/subsets/validation/{{dataset}}_pos.csv",
        validation_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_neg.csv",
        test_pos = f"work_folder{pn}/subsets/test/{{dataset}}_pos.csv",
        test_neg = f"work_folder{pn}/subsets/test/{{dataset}}_neg.csv",
    script:
        "scripts/define_train_test.py"
