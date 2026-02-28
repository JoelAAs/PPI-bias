import networkx as nx
import numpy as np
import argparse
import pandas as pd



def get_negative_data_undir(G_pos):
    nodes = list(G_pos.nodes())
    n_edges = len(G_pos.edges())

    degrees = dict(G_pos.degree())
    probs = np.array([degrees[n] for n in nodes], dtype=float)
    probs /= probs.sum()

    pos_edges = {tuple(sorted((u, v))) for u, v in G_pos.edges()}
    chosen_edges = set()
    

    n_random_edges = 0
    chosen_edges = []
    while n_random_edges < n_edges:
        u = np.random.choice(nodes, p=probs)
        v = np.random.choice(nodes, p=probs)
        
        if u == v:
            continue
        
        neg_edge = tuple(sorted((u,v)))
        
        if neg_edge not in pos_edges and neg_edge not in chosen_edges:
            chosen_edges.append(neg_edge)
            n_random_edges += 1
    
     
    return chosen_edges   ## to dataframe edgelist with node name for each edge


def get_negative_data_dir(G_dir_pos):
    n_edges = len(G_dir_pos.edges())

    prey_degrees = dict(G_dir_pos.out_degree())
    prey_nodes = list(prey_degrees.keys()) 
    prey_probs = np.array([prey_degrees[n] for n in prey_nodes], dtype=float)
    prey_probs /= prey_probs.sum()

    bait_degrees = dict(G_dir_pos.out_degree())
    bait_nodes = list(bait_degrees.keys())
    bait_probs = np.array([bait_degrees[n] for n in bait_nodes], dtype=float)
    bait_probs /= bait_probs.sum()

    pos_edges = {(u, v) for u, v in G_dir_pos.edges()}
    chosen_edges = set()
    

    n_random_edges = 0
    chosen_edges = []
    while n_random_edges < n_edges:
        bait_u = np.random.choice(bait_nodes, p=bait_probs)
        prey_v = np.random.choice(prey_nodes, p=prey_probs)
        
        if bait_u == prey_v:
            continue
        
        neg_edge = (bait_u,prey_v)
        
        if neg_edge not in pos_edges and neg_edge not in chosen_edges:
            chosen_edges.append(neg_edge)
            n_random_edges += 1
    
    return chosen_edges   ## to dataframe edgelist with node name for each edge




def write_edges(edges, flip, output):
    df_neg = pd.DataFrame(edges, columns=["gene_name_bait", "gene_name_prey"])
    if flip:
        flipped_df = df_neg.copy()
        flipped_df[["gene_name_bait", "gene_name_prey"]] = (
            flipped_df[["gene_name_prey", "gene_name_bait"]]
        )

        df_neg = (
            pd.concat([df_neg, flipped_df], ignore_index=True)
        )
    df_neg.to_parquet(output, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Balance positive and negative PPI datasets by matching degree distributions.")
    parser.add_argument("--positive_data", type=str, required=True, help="Path to the positive PPI dataset (CSV format).")
    parser.add_argument("--negative_data", type=str, required=True, help="Path to the negative PPI output.")
    parser.add_argument("--network_type", type=str, required=True, help="Networktype either directed or undirected")
    
    args = parser.parse_args()
    df_pos = pd.read_parquet(args.positive_data)
    G_pos = get_positive_graph(df_pos)
    if args.network_type == "undirectional":
        G_pos = nx.from_pandas_edgelist(df_pos, source="gene_name_bait", target="gene_name_prey", create_using=nx.Graph())
        edges = get_negative_data_undir(G_pos)
        write_edges(edges, True, args.negative_data)
    elif args.networktype == "directional":
        G_pos = nx.from_pandas_edgelist(df_pos, source="gene_name_bait", target="gene_name_prey", create_using=nx.DiGraph())
        edges = get_negative_data_dir(G_pos)
        write_edges(edges, False, args.negative_data)
    else:
        raise ValueError("Wrong")
    