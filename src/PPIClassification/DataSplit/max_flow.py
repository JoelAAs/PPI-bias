import networkx as nx
from fractions import Fraction
import pandas as pd
import argparse
import numpy as np
from scipy.stats import spearmanr

def get_degrees(G):
    degree_bait = dict(G.out_degree())
    degree_prey = dict(G.in_degree())
    
    return degree_bait, degree_prey


def get_scaled_targets(graph, scale):
    bait_target, prey_target = get_degrees(graph)

    bait_target = {gene: round(degree * scale) for gene, degree in bait_target.items()}
    prey_target = {gene: round(degree * scale) for gene, degree in prey_target.items()}

    return bait_target, prey_target


def build_flow_graph(graph, target_bait, target_prey):
    F = nx.DiGraph()

    for node in graph.nodes():
        if node in target_bait:
            F.add_edge("source", ("bait", node), capacity=target_bait.get[node])
        if node in target_prey:
            F.add_edge(("prey", node), "sink", capacity=target_prey[node])

    for bait, prey in graph.edges():
        F.add_edge(("bait", bait), ("prey", prey), capacity=1)

    return F


def get_selected_negative_graph(flow_dict, negative_graph):
    S = nx.DiGraph()
    S.add_nodes_from(negative_graph.nodes())

    for u, v in negative_graph.edges():
        if flow_dict.get(("out", u), {}).get(("in", v), 0) == 1:
            S.add_edge(u, v)
    return S



def get_degree_divergence(targetG, otherG):
    bait_target, prey_target = get_degrees(targetG)
    bait_selected, prey_selected = get_degrees(otherG)
    all_nodes = set(targetG.nodes()) | set(otherG.nodes())
    
    divergence_bait = 0
    divergence_prey = 0
    
    for node in all_nodes:
        divergence_bait += np.abs(bait_target.get(node, 0) - bait_selected.get(node,0))
        divergence_prey += np.abs(prey_target.get(node, 0) - prey_selected.get(node,0))
    
    return divergence_bait, divergence_prey

def degree_spearman_correlation(targetG, otherG):
    bait_target, prey_target = get_degrees(targetG)
    bait_selected, prey_selected = get_degrees(otherG)
    
    targets = []
    selected = []
    all_nodes = set(targetG.nodes()) | set(otherG.nodes())
    for node in all_nodes:
        if node in bait_selected or node in bait_target:
           targets.append(bait_target.get(node, 0))
           selected.append(bait_selected.get(node, 0))
        
        if node in prey_selected or node in prey_target:
           targets.append(prey_target.get(node, 0))
           selected.append(prey_selected.get(node, 0))
        
        rho, _ = spearmanr(targets, selected)
        return rho



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--max_flow_positive", required=True, help="Path to output csv file")
    parser.add_argument("--max_flow_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balance_file", required=True, help="")
    
    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data
    max_flow_positive = args.max_flow_positive
    max_flow_negative = args.max_flow_negative

    positive_bait_prey_df = pd.read_parquet(positive_data)
    negative_bait_prey_df = pd.read_parquet(negative_data)

    positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
    negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

    positive_bait_prey_df.columns = ["bait", "prey"]
    negative_bait_prey_df.columns = ["bait", "prey"]

    shared_baits = set(negative_bait_prey_df["bait"]) & set(positive_bait_prey_df["bait"])
    shared_prey = set(negative_bait_prey_df["prey"]) & set(positive_bait_prey_df["prey"])

    negative_bait_prey_df = negative_bait_prey_df[
        (negative_bait_prey_df["bait"].isin(shared_baits)) &
        (negative_bait_prey_df["prey"].isin(shared_prey))
    ]

    positive_bait_prey_df = positive_bait_prey_df[
        (positive_bait_prey_df["bait"].isin(shared_baits)) &
        (positive_bait_prey_df["prey"].isin(shared_prey))
    ]


    positive_diG = nx.from_pandas_edgelist(
        positive_bait_prey_df, "bait", "prey", create_using=nx.DiGraph()
    )
    negative_diG = nx.from_pandas_edgelist(
        negative_bait_prey_df, "bait", "prey", create_using=nx.DiGraph()
    )
    success = False

    current_min_score = 1000
    n_pos_edges = positive_bait_prey_df.shape[0]
    best_sample_balance = 5
    best_divergence_balance = 1
    best_spearman = 0
    current_best_negative = None
    with open(args.balance_file, "w") as w:
        w.write("positive_edges\tnegative_edges\tscale\tspearman_degree\tdivergance_bait\tdivergrence_prey\n")
        for scale in np.linespace(1, 2, 0.1):
            target_in, target_out = get_scaled_targets(positive_diG, scale)
            print(f"Trying a subset with scaling: {scale}")
            F = build_flow_graph(negative_diG, target_in, target_out)
            flow_value, flow_dict = nx.maximum_flow(F, "source", "sink")
            
            percent_output = round(flow_value / sum(target_in.values()).numerator * 100)

            selected_negative_G = get_selected_negative_graph(flow_dict, negative_diG)
            spearman = degree_spearman_correlation(positive_diG, selected_negative_G)
            div_bait, div_prey = get_degree_divergence(positive_diG, selected_negative_G)

            n_negative_edges = len(list(selected_negative_G.edges()))

            score = 1 - np.abs(n_pos_edges/n_negative_edges) + (div_bait+div_prey)/n_pos_edges - spearman
            w.write(f"{n_pos_edges}\t{n_negative_edges}\t{scale}\t{spearman}\t{div_bait}\t{div_prey}\n")
            if score < current_min_score:
                current_best_negative = selected_negative_G

    
    for balanced_network, output_filename in zip([positive_diG, current_best_negative], [max_flow_positive, max_flow_negative]):
        with open(output_filename, "w") as w:
            for bait, prey in balanced_network.edges():
                w.write(f"{bait}\t{prey}\n")
