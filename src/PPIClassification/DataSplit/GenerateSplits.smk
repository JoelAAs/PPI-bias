import pandas as pd
import itertools
import random
import networkx as nx
from scripts.sample_balance_multi_network_functions import filter_on_threshold

rule prepare_classification_source_max_min:
    params:
        source = lambda wildcards: config["classification_sources"][wildcards.source],
        bait_col = f"{config['id_pattern']}_bait",
        prey_col = f"{config['id_pattern']}_prey"
    input:
        table = lambda wildcards: config["classification_sources"][wildcards.source]["path"].format(
            network_type=wildcards.network_type, dataset=wildcards.dataset)
    output:
        prepared = "work_folder/subsets/{source}/input/{dataset}_{network_type}.pq",
        max_positive = "work_folder/subsets/{source}/strict/{dataset}_{network_type}_pos.csv",
        max_negative = "work_folder/subsets/{source}/strict/{dataset}_{network_type}_neg.csv",
    resources:
        mem_gb=64
    log:
        "logs/subsets/{source}/input/{dataset}_{network_type}.log"
    run:
        positive, negative = params.source["positive"], params.source["negative"]
        value_cols = list(set(["n_observed", positive["column"], negative["column"]])) # if posterior
        source_df = pd.read_parquet(input.table, columns=[params.bait_col, params.prey_col] + value_cols)
        source_df = source_df.rename(columns={params.bait_col: "bait", params.prey_col: "prey"})
        source_df.to_parquet(output.prepared, index=False)

        positive_edges = filter_on_threshold(source_df, positive, positive["strict"])
        negative_edges = filter_on_threshold(source_df, negative, negative["strict"])

        positive_edges[["bait", "prey"]].to_csv(output.max_positive, index=False, sep="\t")
        negative_edges[["bait", "prey"]].to_csv(output.max_negative, index=False, sep="\t")


rule generate_test_validation_test:
    input:
        interaction_data = "work_folder/subsets/{source}/strict/{dataset}_{network_type}_pos.csv",
        max_negative = "work_folder/subsets/{source}/strict/{dataset}_{network_type}_neg.csv"
    output:
        validation_pos = "work_folder/subsets/{source}/validation/{dataset}_{network_type}_pos.csv",
        validation_neg = "work_folder/subsets/{source}/validation/{dataset}_{network_type}_neg.csv",
        test_pos = "work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_neg = "work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
    log:
        "logs/subsets/{source}/validation/{dataset}_{network_type}.log"
    script:
        "scripts/define_validation_test.py"


rule generate_threshold_datasets:
    params:
        source = lambda wildcards: config["classification_sources"][wildcards.source]
    input:
        test_set = "work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        validation_set = "work_folder/subsets/{source}/validation/{dataset}_{network_type}_pos.csv",
        full_detection = "work_folder/subsets/{source}/input/{dataset}_{network_type}.pq"
    threads: 10
    resources:
        mem_gb=130
    output:
        balanced_edges_positive = "work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        balanced_edges_negative = "work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_neg.csv"
    log:
        "logs/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_balance.log"
    script:
        "scripts/single_balanced_threshold_datasets.py"



rule generate_negative_sample:
    input:
        balanced_positive_edges = "work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
    output:
        random_negative_edges = "work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-random_neg.csv"
    log:
        "logs/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_random_neg.log"
    resources:
        mem_gb=80
    script:
        "scripts/random_negative.py"
