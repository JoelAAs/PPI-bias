import networkx as nx
from fractions import Fraction
import pandas as pd
import argparse
import numpy as np
from scipy.stats import spearmanr
from graph_tool.all import Graph, openmp_set_num_threads
from graph_tool.flow import boykov_kolmogorov_max_flow as max_flow


def generate_graph(edge_df, node_map):
    g = Graph(directed=True)
    g.add_vertex(len(node_map))


    bait_src = edge_df["bait"].map(node_map).to_numpy()
    tar_prey = edge_df["prey"].map(node_map).to_numpy()
    edges = np.column_stack((bait_src, tar_prey))

    g.add_edge_list(edges)
    
    return g

def get_degree(g):
    # bait, prey degree
    return g.get_out_degrees(g.get_vertices()), g.get_in_degrees(g.get_vertices())


def get_scaled_targets(graph, scale):
    bait_target, prey_target = get_degree(graph)

    bait_target = np.round(bait_target * scale).astype(np.int64)
    prey_target = np.round(prey_target * scale).astype(np.int64)

    return bait_target, prey_target

def get_node_map(all_nodes):
    all_nodes = list(all_nodes)
    return {gene:i for i, gene in enumerate(all_nodes)}, {i:gene for i, gene in enumerate(all_nodes)}


def build_flow_graph_gt(negative_df, target_bait, target_prey, node_map):
    # 0 is bait, 1 is prey in tuples
    g = Graph(directed=True)

    flow_node_map = {(0,i):i for i, _ in enumerate(target_bait)}
    flow_node_map.update({(1,i):(i+ len(target_bait)) for i, _ in enumerate(target_prey)})
    flow_node_map_index = {value:key for key, value in flow_node_map.items()}
    
    g.add_vertex(len(flow_node_map))

    source = g.add_vertex()
    sink = g.add_vertex()


    capacity = g.new_edge_property("int")

    for bait, cap in enumerate(target_bait):
        bait_v = flow_node_map[(0, bait)]
        e = g.add_edge(source, bait_v)
        capacity[e] = int(cap)

    for prey, cap in enumerate(target_prey):
        prey_v = flow_node_map[(1, prey)]
        e = g.add_edge(prey_v, sink)
        capacity[e] = int(cap)

    src_baits = negative_df["bait"].map(node_map).to_numpy()
    tar_prey = negative_df["prey"].map(node_map).to_numpy()

    bait_ids = [flow_node_map[(0, i)] for i in src_baits]
    prey_ids = [flow_node_map[(1, i)] for i in tar_prey]

    edges = np.column_stack((bait_ids, prey_ids))
    g.add_edge_list(edges)
    for e in list(g.edges())[-len(edges):]:
        capacity[e] = 1
    
    return g, capacity, source, sink, flow_node_map_index


def extract_selected_edges(flow_g, capacity, residual, flow_node_map_index, node_map_index, source, sink):
    selected_edges = []

    for e in flow_g.edges():
        u = e.source()
        v = e.target()

        if u == source or v == sink:
            continue

        flow = capacity[e] - residual[e] # if not in residual then edge is chosen

        if flow == 1:
            bait = node_map_index[flow_node_map_index[int(u)][1]]
            prey = node_map_index[flow_node_map_index[int(v)][1]]
            selected_edges.append((bait, prey))

    return pd.DataFrame(selected_edges, columns=["bait", "prey"])



def get_degree_divergence(targetG, otherG, node_map):
    bait_target, prey_target = get_degree(targetG)
    bait_selected, prey_selected = get_degree(otherG)
    
    divergence_bait = 0
    divergence_prey = 0
    
    for node in node_map.values():
        divergence_bait += np.abs(bait_target[node] - bait_selected[node])
        divergence_prey += np.abs(prey_target[node] - prey_selected[node])
    
    return divergence_bait, divergence_prey

def degree_spearman_correlation(targetG, otherG):
    bait_target, prey_target = get_degree(targetG)
    bait_selected, prey_selected = get_degree(otherG)

    targets = np.concatenate([bait_target, prey_target])
    selected = np.concatenate([bait_selected, prey_selected])

    rho, _ = spearmanr(targets, selected)
    return rho



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--max_flow_positive", required=True, help="Path to output csv file")
    parser.add_argument("--max_flow_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balance_file", required=True, help="")
    parser.add_argument("--threads", type=int, default=10)
    

    args = parser.parse_args()
    openmp_set_num_threads(args.threads)
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



    node_map, node_map_idx = get_node_map(shared_baits | shared_prey)

    pos_diG = generate_graph(positive_bait_prey_df, node_map)

    current_min_score = 1000
    n_pos_edges = positive_bait_prey_df.shape[0]
    best_sample_balance = 5
    best_divergence_balance = 1
    best_spearman = 0 
    current_best_negative = None


    with open(args.balance_file, "w") as w:
        w.write("positive_edges\tnegative_edges\tscale\tspearman_degree\tdivergence_bait\tdivergrence_prey\n")
        for scale in np.linspace(1, 2, 10):

            target_bait, target_prey = get_scaled_targets(pos_diG, scale)

            g, capacity, source, sink, flow_node_map_idx = build_flow_graph_gt(
                negative_bait_prey_df,
                target_bait,
                target_prey,
                node_map
            )
            residual = max_flow(g, source, sink, capacity)

            flow_value = sum(
                capacity[e] - residual[e]
                for e in source.out_edges()
            )

            selected_negative_edges = extract_selected_edges(
                g,
                capacity,
                residual,
                flow_node_map_idx,
                node_map_idx,
                source, 
                sink
            )

            selected_negative_g = generate_graph(selected_negative_edges, node_map)

            
            percent_output = round(flow_value / sum(target_bait) * 100)

            spearman = degree_spearman_correlation(pos_diG, selected_negative_g)
            div_bait, div_prey = get_degree_divergence(pos_diG, selected_negative_g, node_map)

            n_negative_edges = selected_negative_edges.shape[0]

            score = 1 - np.abs(n_pos_edges/n_negative_edges) + (div_bait+div_prey)/n_pos_edges + 1 - spearman
            w.write(f"{n_pos_edges}\t{n_negative_edges}\t{scale}\t{spearman}\t{div_bait}\t{div_prey}\n")
            if score < current_min_score:
                current_min_score = score
                current_best_negative = selected_negative_g

    
    for balanced_network, output_filename in zip([pos_diG, current_best_negative], [max_flow_positive, max_flow_negative]):
        with open(output_filename, "w") as w:
            for bait, prey in balanced_network.edges():
                w.write(f"{node_map_idx[bait]}\t{node_map_idx[prey]}\n")
