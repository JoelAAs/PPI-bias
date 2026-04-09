from degree_balancing_functions import *
import pandas as pd
import networkx as nx


def define_validation_test(G_pos, G_neg, max_iterations):
    overlapping_G, non_overlapping_nodes = remove_all_nonovelapping_nodes(
        G_pos.copy(), G_neg.copy()
    )
    cG_pos, cG_neg = overlapping_G

    remaining_nodes = list(set(cG_pos.nodes()) | set(cG_neg.nodes()))
    n_validate_nodes = round(len(remaining_nodes) * 0.9)
    current_iteration = 0

    best_score = np.inf
    max_tries = 5
    tries = 0
    prev_size_balance = None

    while current_iteration < max_iterations:
        current_iteration += 1
        validation_nodes = set()
        test_nodes = non_overlapping_nodes.copy()
        classes = np.zeros_like(remaining_nodes, dtype=int)
        classes[:n_validate_nodes] = 1
        np.random.shuffle(classes)
        for c, node in zip(classes, remaining_nodes):
            if c == 0:
                test_nodes.add(node)
            else:
                validation_nodes.add(node)

        rG_validate_pos, rG_validate_neg = cG_pos.subgraph(
            validation_nodes
        ), cG_neg.subgraph(validation_nodes)
        G_validate_pos, G_validate_neg, discarded_nodes = back_and_forth_max_flow(
            rG_validate_pos, rG_validate_neg
        )
        test_nodes |= discarded_nodes

        G_test_pos, G_test_neg = G_pos.subgraph(test_nodes), G_neg.subgraph(test_nodes)
        # validate_balance = np.abs(G_validate_pos.number_of_edges() - G_validate_neg.number_of_edges()) always gonna be fine
        size_balance = G_validate_pos.number_of_edges() - G_test_pos.number_of_edges()

        step_size = size_balance / (
            G_validate_pos.number_of_edges() + G_test_pos.number_of_edges()
        )
        next_node_count = n_validate_nodes - round(step_size * len(classes))
        if next_node_count < 1:
            n_validate_nodes = 0
        elif next_node_count > len(remaining_nodes):
            n_validate_nodes = len(remaining_nodes)
        else:
            n_validate_nodes = next_node_count

        if np.abs(size_balance) < best_score:
            best_test = [G_test_pos.copy(), G_test_neg.copy()]
            best_validation = [G_validate_pos.copy(), G_validate_neg.copy()]
            best_score = np.abs(size_balance)

        if size_balance == prev_size_balance:
            tries += 1
        prev_size_balance = size_balance
        if tries > max_tries:
            break

    return *best_validation, *best_test


def main():
    hci_df = pd.read_csv(snakemake.input.interaction_data, sep="\t")[
        ["gene_name_bait", "gene_name_prey"]
    ]
    hcni_df = pd.read_csv(snakemake.input.max_negative, sep="\t")[
        ["gene_name_bait", "gene_name_prey"]
    ]

    G_pos = nx.from_pandas_edgelist(
        hci_df, "gene_name_bait", "gene_name_prey", create_using=nx.DiGraph
    )
    G_neg = nx.from_pandas_edgelist(
        hcni_df, "gene_name_bait", "gene_name_prey", create_using=nx.DiGraph
    )

    G_validation_pos, G_validation_neg, G_test_pos, G_test_neg = define_validation_test(
        G_pos, G_neg, 500
    )

    nx.write_edgelist(G_validation_pos, snakemake.output.validation_pos, delimiter="\t", data=False)
    nx.write_edgelist(G_validation_neg, snakemake.output.validation_neg, delimiter="\t", data=False)
    nx.write_edgelist(G_test_pos, snakemake.output.test_pos, delimiter="\t", data=False)
    nx.write_edgelist(G_test_neg, snakemake.output.test_neg, delimiter="\t", data=False)
