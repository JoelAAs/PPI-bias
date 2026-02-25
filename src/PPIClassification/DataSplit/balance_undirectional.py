import networkx as nx
import pandas as pd
import numpy as np
import argparse
np.random.seed(1234)


def get_network_from_csv(edege_df):
    G = nx.from_pandas_edgelist(edege_df, source="gene_name_bait", target="gene_name_prey", create_using=nx.Graph())
    return G


def balance_split(df_pos, df_neg, output_pos, output_neg):
    edge_df_pos = pd.read_parquet(df_pos)
    edge_df_neg = pd.read_parquet(df_neg)

    pos_genes = set(edge_df_pos["gene_name_bait"]) | set(edge_df_pos["gene_name_prey"])
    negative_genes = set(edge_df_neg["gene_name_bait"]) | set(edge_df_neg["gene_name_prey"])

    non_shared_genes = (pos_genes | negative_genes) - (pos_genes & negative_genes)

    edge_df_pos = edge_df_pos[(~edge_df_pos["gene_name_bait"].isin(non_shared_genes)) | (~edge_df_pos["gene_name_prey"].isin(non_shared_genes))]
    edge_df_neg = edge_df_neg[(~edge_df_neg["gene_name_bait"].isin(non_shared_genes)) | (~edge_df_neg["gene_name_prey"].isin(non_shared_genes))]

    G_pos = get_network_from_csv(edge_df_pos)
    G_neg = get_network_from_csv(edge_df_neg)
    
    selected_edges = []
    expected_positive_degree = dict(G_pos.degree())
    negative_edges = list(G_neg.edges())
    np.random.shuffle(negative_edges)

    for (u,v) in negative_edges:
        capacity_a = expected_positive_degree.get(u, False)
        capacity_b = expected_positive_degree.get(v, False)
        
        if capacity_a > 0 and capacity_b > 0:
            selected_edges.append((u, v))
            expected_positive_degree[u] -= 1
            expected_positive_degree[v] -= 1
        else:
            
            if capacity_a > 0 and capacity_b == False:
                expected_positive_degree[u] -= 1
                selected_edges.append((u, v))
            elif capacity_b > 0 and capacity_a == False:
                expected_positive_degree[v] -= 1
                selected_edges.append((u, v))

    selected_edges_df = pd.DataFrame(selected_edges, columns=["gene_name_bait", "gene_name_prey"])
    selected_edges_df.columns = ["bait", "prey"]
    selected_edges_df.to_csv(output_neg, sep="\t", index=False)
    edge_df_pos = edge_df_pos[["gene_name_bait", "gene_name_prey"]]
    edge_df_pos.columns = ["bait", "prey"]
    edge_df_pos.to_csv(output_pos, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Balance positive and negative PPI datasets by matching degree distributions.")
    parser.add_argument("--positive_data", type=str, required=True, help="Path to the positive PPI dataset (CSV format).")
    parser.add_argument("--negative_data", type=str, required=True, help="Path to the negative PPI dataset (CSV format).")
    parser.add_argument("--output_positive", type=str, required=True, help="Path to save the balanced positive PPI dataset (CSV format).")
    parser.add_argument("--output_negative", type=str, required=True, help="Path to save the balanced negative PPI dataset (CSV format).")

    args = parser.parse_args()
    balance_split(args.positive_data, args.negative_data, args.output_positive, args.output_negative)