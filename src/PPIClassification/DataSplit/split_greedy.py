import numpy as np
import pandas as pd
import networkx as nx

####################################################################
#
# Split edges while keeping the balance and degree somewhat intact
#
####################################################################


def get_worst_node(G_pos, G_neg):
    bait_degree = [dict(G_pos.out_degree()), dict(G_neg.out_degree())]
    prey_degree = [dict(G_pos.in_degree()), dict(G_neg.in_degree())]

    worst_node = None
    worst_score = 0

    def _get_delta(b_degree, p_degree):
        score = 0
        for diff_set in [b_degree, p_degree]:
            min_deg = min(diff_set)
            delta = np.abs(diff_set[0] - diff_set[1])
            if sum(diff_set) == 0:
                continue
            elif min_deg == 0:
                score += 100000  # Large number to always pick those without negative/positive representation
            else:
                score += delta / min_deg

        return score

    for node in G_pos:
        b_n_degree = [bait_degree[0][node], bait_degree[1][node]]
        p_n_degree = [prey_degree[0][node], prey_degree[1][node]]

        score = _get_delta(b_n_degree, p_n_degree)
        if score > worst_score:
            worst_score = score
            worst_node = node

    return worst_node


def get_bait_prey_div(G_pos, G_neg):
    bait_degree = [dict(G_pos.out_degree()), dict(G_neg.out_degree())]
    prey_degree = [dict(G_pos.in_degree()), dict(G_neg.in_degree())]

    df_delta = pd.DataFrame({"node": G_pos.nodes()})
    df_delta["pos_bait"] = df_delta["node"].apply(lambda x: bait_degree[0][x])
    df_delta["neg_bait"] = df_delta["node"].apply(lambda x: bait_degree[1][x])
    df_delta["pos_prey"] = df_delta["node"].apply(lambda x: prey_degree[0][x])
    df_delta["neg_prey"] = df_delta["node"].apply(lambda x: prey_degree[1][x])

    df_delta["delta_bait"] = (df_delta["pos_bait"] - df_delta["neg_bait"]).abs()
    df_delta["delta_prey"] = (df_delta["pos_prey"] - df_delta["neg_prey"]).abs()

    return df_delta["delta_bait"].sum(), df_delta["delta_prey"].sum()


def report_balance(G_pos, G_neg):
    bait_div, prey_div = get_bait_prey_div(G_pos, G_neg)
    msg = (
        f"------------------------------------------\n"
        f"Sample balance {G_pos.number_of_edges()}/{G_neg.number_of_edges()} ({round(G_pos.number_of_edges()/ G_neg.number_of_edges()*100)})\n"
        f"Bait div: {bait_div} ({bait_div/ G_pos.number_of_edges()})\t Prey div {prey_div} ({prey_div/ G_pos.number_of_edges()})\n"
    )

    print(msg, flush=True)


def remove_nodes_until_edge_count(G_pos, G_neg, fraction_to_pick):

    G_pos_train = G_pos.copy()
    G_neg_train = G_neg.copy()

    removed_nodes = []
    i = 1
    current_ratio = 1
    old_ratio = current_ratio
    while current_ratio > fraction_to_pick:
        i += 1
        w_node = get_worst_node(G_pos_train, G_neg_train)
        removed_nodes.append(w_node)
        G_pos_train.remove_node(w_node)
        G_neg_train.remove_node(w_node)

        G_pos_remain = G_pos.subgraph(removed_nodes)
        G_neg_remain = G_neg.subgraph(removed_nodes)

        total_remaining_edges = (
            G_pos_remain.number_of_edges() + G_pos_train.number_of_edges()
        )
        current_ratio = G_pos_train.number_of_edges() / total_remaining_edges
        if old_ratio - current_ratio > 0.01:
            print(
                f"Current progress: {round(current_ratio*100)} %; Requested ratio {round(fraction_to_pick*100)} %",
                flush=True,
            )
            old_ratio = current_ratio

    print("Original set")
    report_balance(G_pos, G_neg)
    print("Selected set")
    report_balance(G_pos_train, G_neg_train)
    print("Remaining edges")
    report_balance(G_pos_remain, G_neg_remain)
    return [G_pos_train, G_neg_train], [G_pos_remain, G_neg_remain]


####################################################################
#
# Balance with bouncy maxflow
#
####################################################################


def get_graph_from_selected_edges(flow_G, idx_to_gene):
    S = nx.DiGraph()

    for b, ps in flow_G.items():
        if b in idx_to_gene:
            for p, flow in ps.items():
                if p in idx_to_gene and flow == 1:
                    S.add_edge(idx_to_gene[b], idx_to_gene[p])

    return S


def get_targets(G, node_idx):
    target_bait = {
        node_idx[gene]: degree_count
        for gene, degree_count in dict(G.out_degree()).items()
    }
    target_prey = {
        node_idx[gene]: degree_count
        for gene, degree_count in dict(G.in_degree()).items()
    }
    return target_bait, target_prey


def build_flow_graph(target_bait, target_prey, edge_G, node_idx):
    F = nx.DiGraph()
    bp_index = {(0, b_idx): i for i, b_idx in enumerate(target_bait.keys())}
    bp_index.update(
        {(1, p_idx): i + len(target_bait) for i, p_idx in enumerate(target_prey.keys())}
    )
    source = len(bp_index)
    sink = len(bp_index) + 1
    for node, degree_val in target_bait.items():
        bait_node = bp_index[(0, node)]
        F.add_edge(source, bait_node, capacity=degree_val)

    for node, degree_val in target_prey.items():
        prey_node = bp_index[(1, node)]
        F.add_edge(prey_node, sink, capacity=degree_val)

    for b, p in edge_G.edges():
        b_idx = node_idx[b]
        p_idx = node_idx[p]
        if (0, b_idx) in bp_index and (1, p_idx) in bp_index:
            bait_node = bp_index[(0, b_idx)]
            prey_node = bp_index[(1, p_idx)]
            F.add_edge(bait_node, prey_node, capacity=1)

    return F, source, sink, bp_index


def remove_all_nonovelapping_nodes(G_pos, G_neg):
    any_nonoverlapping = True
    nodes_to_drop = set()

    while any_nonoverlapping:
        for node in nodes_to_drop:
            for G in [G_pos, G_neg]:
                try:
                    G.remove_node(node)
                except nx.NetworkXError:
                    continue
        if G_pos.number_of_edges() == 0 or G_neg.number_of_edges() == 0:
            raise nx.NetworkXError("Graph is empty!")
        
        pos_bait = {u for u, _ in G_pos.edges()}
        neg_bait = {u for u, _ in G_neg.edges()}
        pos_prey = {v for _, v in G_pos.edges()}
        neg_prey = {v for _, v in G_neg.edges()}

        non_overlapping_baits = pos_bait ^ neg_bait
        non_overlapping_prey = pos_prey ^ neg_prey

        nodes_to_drop = non_overlapping_baits | non_overlapping_prey

        if len(nodes_to_drop) == 0:
            any_nonoverlapping = False

    return G_pos, G_neg


def back_and_forth_max_flow(G_pos, G_neg):
    i = 0
    target_source = {"pos": G_pos, "neg": G_neg}

    node_idx = {
        gene: i for i, gene in enumerate(set(G_pos.nodes()) | set(G_neg.nodes()))
    }
    node_idx_gene = {i: gene for gene, i in node_idx.items()}
    percent_flow = 0
    keys = list(target_source.keys())
    while percent_flow < 0.9:
        target_key = keys[i % 2]
        source_key = keys[(i + 1) % 2]
        target_G = target_source[target_key]
        edge_G = target_source[source_key]

        target_bait, target_prey = get_targets(target_G, node_idx)

        mf_F, source, sink, flow_idx = build_flow_graph(
            target_bait, target_prey, edge_G, node_idx
        )
        flow_value, flow_dict = nx.maximum_flow(mf_F, source, sink)
        flow_idx_gene = {i: node_idx_gene[idx[1]] for idx, i in flow_idx.items()}

        percent_flow = flow_value / sum(target_bait.values())
        print(f"Current_flow is {percent_flow}")

        selected_G = get_graph_from_selected_edges(flow_dict, flow_idx_gene)
        balanced_source = remove_all_nonovelapping_nodes(target_G, selected_G)

        target_source = {target_key: balanced_source[0], source_key: balanced_source[1]}
        i += 1

    return target_source["pos"], target_source["neg"]


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

    train_graphs, remaining_graphs = remove_nodes_until_edge_count(G_pos, G_neg, 0.7)
    G_pos_train, G_neg_train = back_and_forth_max_flow(train_graphs[0], train_graphs[1])

    validation_graphs, test_graphs = remove_nodes_until_edge_count(
        remaining_graphs[0], remaining_graphs[1], 0.5
    )
    G_pos_validation, G_neg_validation = back_and_forth_max_flow(
        validation_graphs[0], validation_graphs[1]
    )

    G_pos_test, G_neg_test = test_graphs

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

    nx.write_edgelist(
        G_pos_test, snakemake.output.test_pos, delimiter="\t", data=False
    )
    nx.write_edgelist(
        G_neg_test, snakemake.output.test_neg, delimiter="\t", data=False
    )


# Off we go
main()
