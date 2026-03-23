import pandas as pd
import itertools
import random

def get_gene_partition(wc):
    if wc["partition_name"] == "sequencesimilarity":
        return f"work_folder{pn}/subsets/partitions/sequencesimilarity_gene_name.txt"
    elif wc["partition_name"] == "maxpos":
        return f"work_folder{pn}/subsets/partitions/{wc['dataset']}_{wc["network_type"]}_limit_{wc['pos_limit']}_gene_name.txt"
    elif wc["partition_name"] == "pre-balanced":
        return f"work_folder{pn}/subsets/partitions/{wc['dataset']}__directional_limit_{wc['neg_limit']}_poslim_{wc['pos_limit']}_gene_name.txt"
    else:
        raise ValueError(f"unknown partition name = {wc['partition_name']}")


def get_number_of_edges(G, vertices):
    return G.subgraph(vertices).number_of_edges()

def check_fit_of_partition(G_pos, G_neg, partition_list, split_fractions):
    negative_edges = np.array([get_number_of_edges(G_neg, p) for p in partition_list])
    positive_edges = np.array([get_number_of_edges(G_pos, p) for p in partition_list])

    fraction_positive_partition = positive_edges / sum(positive_edges)
    fraction_negative_partition = negative_edges / sum(negative_edges)
    delta_positive_fraction = np.abs(fraction_positive_partition - split_fractions)
    delta_negative_fraction = np.abs(fraction_negative_partition - split_fractions)
    
    fraction_imbalance = sum(delta_positive_fraction+delta_negative_fraction)

    sqr_mean_imbalance = ((negative_edges - positive_edges)/(negative_edges + positive_edges)**2)
    
    total_edges = G_pos.number_of_edges() + G_neg.number_of_edges()
    percent_edges_discarded = 1 - sum(positive_edges + negative_edges) / total_edges

    fitness_score = fraction_imbalance + percent_edges_discarded + sqr_mean_imbalance

    return fitness_score

def get_vertex_list(vertex_blocks, partition_list):
    vertices = {}
    for block in vertex_blocks:
        vertices |= set(partition_list[block])
    return vertices

def get_best_partition(G_pos, G_neg, partition_list, split_fractions, n_shuffles=1000):
    global_best_fit = np.inf
    global_best_set = None
    for _ in range(n_shuffles):
        random.shuffle(partition_list)
        current_best_set = [set() for _ in range(3)]
        best_set = [set() for _ in range(3)]

        for partition in partition_list:
            best_fit = np.inf
            for i in range(3):
                test_set = [s.copy() for s in best_set]
                test_set[i] |= partition
                fit = check_fit_of_partition(G_pos, G_neg, sets, split_fractions)
                if fit < best_fit:
                    best_fit = fit
                    current_best_set = test_set
            best_set = current_best_set
        if best_fit < global_best_fit:
            global_best_fit = best_fit
            global_best_set = best_set

    return global_best_set, global_best_fit


rule define_partitions:
    params:
        split_fractions = [0.7, 0.15, 0.15]
    input:
        gene_partitions=get_gene_partition,
        set_pos=f"work_folder{pn}/subsets/maxflow/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/maxflow/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    output:
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_pos.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_pos.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_neg.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_neg.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_neg.csv"
    run:
        gene_partitions = pd.read_csv(input.gene_partitions, sep="\t")
        pos_df = pd.read_csv(input.set_pos, sep="\t", header=None)
        neg_df = pd.read_csv(input.set_neg, sep="\t", header=None)
        G_pos = nx.from_pandas_edgelist(pos_df, 0, 1)
        G_neg = nx.from_pandas_edgelist(neg_df, 0, 1)
        
        partition_list = (
            gene_partitions
                .groupby("partition")["gene_name"]
                .apply(list)
                .sort_index()
                .tolist()
        )

        best_sets, score = get_best_partition(G_pos, G_neg, partition_list, params.split_fractions)

        # WIP write partition list and get subset of pos and neg edges that fit the partition


