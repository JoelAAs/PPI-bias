import numpy as np
import pandas as pd
from sample_balance_multi_network_functions import drop_exclusive_nodes, degree_balace_edges, sanity_check


def define_test_set(pos_edges, neg_edges, directed, min_flow=0.95):
    original_pos = pos_edges.copy()
    original_neg = neg_edges.copy()

    pos_edges, neg_edges, non_overlapping_nodes = drop_exclusive_nodes(
        pos_edges, neg_edges, directed=directed
    )

    G_test_pos, G_test_neg, discarded_nodes = degree_balace_edges(
        pos_edges, neg_edges, min_flow=min_flow, directed=directed
        )

    return G_test_pos, G_test_neg


def write_edgelist(df, output_file):
    df.to_csv(output_file, sep="\t", index=False)


def main():
    hci_df = pd.read_csv(snakemake.input.interaction_data, sep="\t")
    hrni_df = pd.read_csv(snakemake.input.max_negative, sep="\t")

    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")
    
    G_test_pos, G_test_neg = define_test_set(
        hci_df, hrni_df, directed=directed
    )

    sanity_check(
        [G_test_pos],
        [G_test_neg],
        [hci_df],
        [hrni_df], directed=directed) # validation set
                
    write_edgelist(G_test_pos, snakemake.output.test_pos)
    write_edgelist(G_test_neg, snakemake.output.test_neg)


if __name__ == "__main__":
    main()


