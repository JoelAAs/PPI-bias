from degree_balaning_functions import *


def main():
    pos_df = pd.read_csv(snakemake.input.set_pos, sep="\t", header=None)
    neg_df = pd.read_csv(snakemake.input.set_neg, sep="\t", header=None)

    any_nonoverlapping = True
    while any_nonoverlapping:
        shared_bait = set(pos_df[0]) & set(neg_df[0])
        shared_prey = set(pos_df[1]) & set(neg_df[1])

        pos_df = pos_df[(pos_df[0].isin(shared_bait) & pos_df[1].isin(shared_prey))]
        neg_df = neg_df[(neg_df[0].isin(shared_bait) & neg_df[1].isin(shared_prey))]
        non_overlapping = len(
            (set(pos_df[0]) | set(neg_df[0])) - (set(pos_df[0]) & set(neg_df[0]))
        ) + len((set(pos_df[1]) | set(neg_df[1])) - (set(pos_df[1]) & set(neg_df[1])))
        if non_overlapping == 0:
            any_nonoverlapping = False

    G_pos = nx.from_pandas_edgelist(pos_df, 0, 1, create_using=nx.DiGraph)
    G_neg = nx.from_pandas_edgelist(neg_df, 0, 1, create_using=nx.DiGraph)

    train_graphs, discarded_nodes_ec_train = remove_nodes_until_edge_count(
        G_pos, G_neg, 0.7
    )
    G_pos_train, G_neg_train, discarded_nodes_mf_train = back_and_forth_max_flow(
        train_graphs[0], train_graphs[1]
    )

    discarded_nodes_train = discarded_nodes_ec_train | discarded_nodes_mf_train
    remaining_graphs = [G.subgraph(discarded_nodes_train) for G in [G_pos, G_neg]]

    G_pos_validation, G_neg_validation, G_pos_test, G_neg_test = edge_balance_partition(
        remaining_graphs[0], remaining_graphs[1], max_iter=10000
    )

    nx.write_edgelist(
        G_pos_train, snakemake.output.train_pos, delimiter="\t", data=False
    )
    nx.write_edgelist(
        G_neg_train, snakemake.output.train_neg, delimiter="\t", data=False
    )

    nx.write_edgelist(
        G_pos_validation, snakemake.output.validation_pos, delimiter="\t", data=False
    )
    nx.write_edgelist(
        G_neg_validation, snakemake.output.validation_neg, delimiter="\t", data=False
    )

    nx.write_edgelist(G_pos_test, snakemake.output.test_pos, delimiter="\t", data=False)
    nx.write_edgelist(G_neg_test, snakemake.output.test_neg, delimiter="\t", data=False)


# Off we go
main()
