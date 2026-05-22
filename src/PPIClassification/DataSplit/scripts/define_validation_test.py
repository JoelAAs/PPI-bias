import numpy as np
import pandas as pd
from sample_balance_multi_network_functions import drop_exclusive_nodes, degree_balace_edges, sanity_check


def define_validation_test(pos_edges, neg_edges, max_iterations, directed, min_flow=0.95):
    original_pos = pos_edges.copy()
    original_neg = neg_edges.copy()

    pos_edges, neg_edges, non_overlapping_nodes = drop_exclusive_nodes(
        pos_edges, neg_edges, directed=directed
    )

    remaining_nodes = list(
        set(pos_edges["bait"]) | set(pos_edges["prey"])
        | set(neg_edges["bait"]) | set(neg_edges["prey"])
    )
    n_validate_nodes = round(len(remaining_nodes) * 0.9)

    best_score = np.inf
    best_test = None
    best_validation = None
    max_tries = 5
    tries = 0
    prev_size_balance = None

    for current_iteration in range(1, max_iterations + 1):
        classes = np.zeros(len(remaining_nodes), dtype=int)
        classes[:n_validate_nodes] = 1

        empty_graph_tries = 0
        success = False
        while(not success and empty_graph_tries < 5):
            try:
                validation_nodes = set()
                test_nodes = non_overlapping_nodes.copy()
                np.random.shuffle(classes)
                for c, node in zip(classes, remaining_nodes):
                    if c == 0:
                        test_nodes.add(node)
                    else:
                        validation_nodes.add(node)

                print(f"{len(validation_nodes)} validation nodes")

                val_pos = pos_edges[
                    pos_edges["bait"].isin(validation_nodes) & pos_edges["prey"].isin(validation_nodes)
                ]
                val_neg = neg_edges[
                    neg_edges["bait"].isin(validation_nodes) & neg_edges["prey"].isin(validation_nodes)
                ]

                G_validate_pos, G_validate_neg, discarded_nodes = degree_balace_edges(
                    val_pos, val_neg, min_flow=min_flow, directed=directed
                )
            
                test_nodes |= discarded_nodes

                test_pos = original_pos[
                    original_pos["bait"].isin(test_nodes) & original_pos["prey"].isin(test_nodes)
                ]
                test_neg = original_neg[
                    original_neg["bait"].isin(test_nodes) & original_neg["prey"].isin(test_nodes)
                ]

                size_balance = G_validate_pos.shape[0] - test_pos.shape[0]
                print(
                    f"iteration {current_iteration}: size balance = {size_balance} "
                    f"(Validation: {G_validate_pos.shape[0]}, Test: {test_pos.shape[0]})"
                )
                success= True
            except ValueError as e:
                print(e)
                empty_graph_tries += 1
            
        if (not success and empty_graph_tries >= 5):
            continue

        total = G_validate_pos.shape[0] + test_pos.shape[0]
        if total > 0:
            step_size = size_balance / total
            next_node_count = round(n_validate_nodes * (1 - step_size))
            n_validate_nodes = max(0, min(len(remaining_nodes), next_node_count))

        if np.abs(size_balance) < best_score:
            best_test = [test_pos, test_neg]
            best_validation = [G_validate_pos, G_validate_neg]
            best_score = np.abs(size_balance)

        if size_balance == prev_size_balance:
            tries += 1
        prev_size_balance = size_balance
        if tries > max_tries:
            break

    return *best_validation, *best_test


def write_edgelist(df, output_file):
    df.to_csv(output_file, sep="\t", index=False)


def main():
    hci_df = pd.read_csv(snakemake.input.interaction_data, sep="\t")[
        ["gene_name_bait", "gene_name_prey"]
    ].rename(columns={"gene_name_bait": "bait", "gene_name_prey": "prey"})
    hcni_df = pd.read_csv(snakemake.input.max_negative, sep="\t")[
        ["gene_name_bait", "gene_name_prey"]
    ].rename(columns={"gene_name_bait": "bait", "gene_name_prey": "prey"})

    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")
    
    G_validation_pos, G_validation_neg, G_test_pos, G_test_neg = define_validation_test(
        hci_df, hcni_df, 20, directed=directed
    )

    sanity_check(
        [G_validation_pos],
        [G_validation_neg],
        [hci_df],
        [hcni_df], directed=directed) # validation set

    try: # test-set doesn't need the same node set
        sanity_check(
            [G_test_pos],
            [G_test_neg],
            [hci_df],
            [hcni_df], directed=directed)
    except AssertionError as e:
        if "Node sets do not match" not in str(e):
            raise e
        
    
    write_edgelist(G_validation_pos, snakemake.output.validation_pos)
    write_edgelist(G_validation_neg, snakemake.output.validation_neg)
    write_edgelist(G_test_pos, snakemake.output.test_pos)
    write_edgelist(G_test_neg, snakemake.output.test_neg)


if __name__ == "__main__":
    main()
