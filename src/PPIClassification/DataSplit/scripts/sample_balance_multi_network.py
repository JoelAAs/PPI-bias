import pandas as pd
import numpy as np
from graph_tool.all import Graph, openmp_set_num_threads
from graph_tool.flow import push_relabel_max_flow as max_flow


def generate_graph(edge_df, node_map):
    g = Graph(directed=True)
    g.add_vertex(len(node_map))

    bait_src = edge_df["bait"].map(node_map).to_numpy()
    tar_prey = edge_df["prey"].map(node_map).to_numpy()
    edges = np.column_stack((bait_src, tar_prey))
    edges = edges[~np.isnan(edges).any(axis=1)]
    edges = edges.astype(int)
    g.add_edge_list(edges)

    return g


def get_degree(g):
    return g.get_out_degrees(g.get_vertices()).astype(np.int64), g.get_in_degrees(
        g.get_vertices()
    ).astype(np.int64)


def get_node_map(all_nodes):
    all_nodes = list(all_nodes)
    return {gene: i for i, gene in enumerate(all_nodes)}, {
        i: gene for i, gene in enumerate(all_nodes)
    }


def build_bait_prey_flow_graph(
    edge_df, target_bait, target_prey, g, capacity, node_map
):
    if edge_df.empty:
        raise ValueError("Edge DataFrame is empty. Cannot build flow graph.")
    if any(edge_df["bait"].isin(node_map) == False) or any(
        edge_df["prey"].isin(node_map) == False
    ):
        raise ValueError("Some baits or preys in edge_df are not in node_map.")
    num_vertices = g.num_vertices()

    flow_node_map = {(0, i): (i + num_vertices) for i, _ in enumerate(target_bait)}
    flow_node_map.update(
        {
            (1, i): (i + num_vertices + len(target_bait))
            for i, _ in enumerate(target_prey)
        }
    )
    flow_node_map_index = {value: key for key, value in flow_node_map.items()}

    g.add_vertex(len(flow_node_map))

    source = g.add_vertex()
    sink = g.add_vertex()

    for bait, cap in enumerate(target_bait):
        bait_v = flow_node_map[(0, bait)]
        e = g.add_edge(source, bait_v)
        capacity[e] = int(cap)

    for prey, cap in enumerate(target_prey):
        prey_v = flow_node_map[(1, prey)]
        e = g.add_edge(prey_v, sink)
        capacity[e] = int(cap)

    src_baits = edge_df["bait"].map(node_map).to_numpy(dtype=int)
    tar_prey = edge_df["prey"].map(node_map).to_numpy(dtype=int)

    num_existing = g.num_edges()
    bait_ids = src_baits + num_vertices
    prey_ids = tar_prey + num_vertices + len(target_bait)

    edges = np.column_stack((bait_ids, prey_ids))
    g.add_edge_list(edges)

    # Set capacity=1 for all newly added edges in bulk
    capacity.a[num_existing : num_existing + len(edges)] = 1

    return g, capacity, source, sink, flow_node_map_index, flow_node_map


def extract_selected_edges(
    flow_g, capacity, residual, flow_node_map_index, node_map, source, sink
):
    g = Graph(directed=True)
    g.add_vertex(len(node_map))

    n_verts = flow_g.num_vertices()
    valid = np.zeros(n_verts, dtype=bool)
    orig_idx = np.zeros(n_verts, dtype=int)
    vids = np.array(list(flow_node_map_index.keys()))
    idxs = np.array([v[1] for v in flow_node_map_index.values()])
    valid[vids] = True
    orig_idx[vids] = idxs

    edges = flow_g.get_edges(
        [flow_g.edge_index]
    )  # (n_edges, 3): [source, target, edge_idx]
    edge_flows = capacity.a[edges[:, 2]] - residual.a[edges[:, 2]]

    mask = (edge_flows == 1) & valid[edges[:, 0]] & valid[edges[:, 1]]
    selected = orig_idx[edges[mask, :2]]

    if len(selected):
        g.add_edge_list(selected)
    return g


def graph_to_edge_df(g, node_map_idx):
    edge_list = []
    for e in g.edges():
        bait_idx = int(e.source())
        prey_idx = int(e.target())
        bait = node_map_idx[bait_idx]
        prey = node_map_idx[prey_idx]
        edge_list.append((bait, prey))
    return pd.DataFrame(edge_list, columns=["bait", "prey"])


def single_alternating_maxflow(positive_edge_df, negative_edge_df, min_flow):
    all_baits = set(negative_edge_df["bait"]) & set(positive_edge_df["bait"])
    all_prey = set(negative_edge_df["prey"]) & set(positive_edge_df["prey"])

    node_map, node_map_idx = get_node_map(all_baits | all_prey)
    g_pos = generate_graph(positive_edge_df, node_map)
    target_bait_pos, target_prey_pos = get_degree(g_pos)

    mf_g_pos = Graph(directed=True)
    capacity_pos = mf_g_pos.new_edge_property("int")
    (
        mf_g_pos,
        capacity_pos,
        source_pos,
        sink_pos,
        flow_node_map_index_pos,
        flow_node_map_pos,
    ) = build_bait_prey_flow_graph(
        edge_df=negative_edge_df,
        target_bait=target_bait_pos,
        target_prey=target_prey_pos,
        g=mf_g_pos,
        capacity=capacity_pos,
        node_map=node_map,
    )

    i = 0
    percent_flow_value = 0
    current_g = [mf_g_pos, None]
    current_capacity = [capacity_pos, None]
    current_source = [source_pos, None]
    current_sink = [sink_pos, None]
    current_flow_map_index = [flow_node_map_index_pos, None]
    current_flow_map = [flow_node_map_pos, None]
    subset_edges = [positive_edge_df, negative_edge_df]

    while percent_flow_value < min_flow:
        residual = max_flow(
            current_g[i % 2],
            current_source[i % 2],
            current_sink[i % 2],
            current_capacity[i % 2],
        )

        capacity_edges = [e for e in current_source[i % 2].out_edges()]
        percent_flow_value = sum(
            current_capacity[i % 2][e] - residual[e] for e in capacity_edges
        ) / sum(current_capacity[i % 2][e] for e in capacity_edges)

        print(f"Flow value {percent_flow_value}")

        g_selected = extract_selected_edges(
            flow_g=current_g[i % 2],
            capacity=current_capacity[i % 2],
            residual=residual,
            flow_node_map_index=current_flow_map_index[i % 2],
            node_map=node_map,
            source=current_source[i % 2],
            sink=current_sink[i % 2],
        )
        new_target_bait, new_target_prey = get_degree(g_selected)
        g_next = Graph(directed=True)
        capacity_next = g_next.new_edge_property("int")
        subset_edges[(i + 1) % 2] = graph_to_edge_df(g_selected, node_map_idx)

        (
            current_g[(i + 1) % 2],
            current_capacity[(i + 1) % 2],
            current_source[(i + 1) % 2],
            current_sink[(i + 1) % 2],
            current_flow_map_index[(i + 1) % 2],
            current_flow_map[(i + 1) % 2],
        ) = build_bait_prey_flow_graph(
            edge_df=subset_edges[i % 2],
            target_bait=new_target_bait,
            target_prey=new_target_prey,
            g=g_next,
            capacity=capacity_next,
            node_map=node_map,
        )
        i += 1

    return subset_edges[0], subset_edges[1], node_map_idx


def build_multi_network_flow_graph(edge_list_a, edge_list_b):
    if len(edge_list_a) != len(edge_list_b):
        raise ValueError("Positive and negative edge lists must have the same length.")
    min_flow = 0.9

    smallest_edge_count = min(
        single_alternating_maxflow(a_edge_df, b_edge_df, min_flow)[0].shape[0]
        for a_edge_df, b_edge_df in zip(edge_list_a, edge_list_b)
    )
    network_data = dict()
    multi_flow_graph = Graph(directed=True)
    mf_capacity = multi_flow_graph.new_edge_property("int")
    super_source = multi_flow_graph.add_vertex()
    super_sink = multi_flow_graph.add_vertex()

    for net_i, (a_edge_df, b_edge_df) in enumerate(zip(edge_list_a, edge_list_b)):
        shared_baits = set(a_edge_df["bait"]) & set(b_edge_df["bait"])
        shared_prey = set(a_edge_df["prey"]) & set(b_edge_df["prey"])

        node_map, node_map_idx = get_node_map(shared_baits | shared_prey)
        g_a = generate_graph(a_edge_df, node_map)
        target_bait, target_prey = get_degree(g_a)
        (
            multi_flow_graph,
            mf_capacity,
            source,
            sink,
            flow_node_map_index,
            flow_node_map,
        ) = build_bait_prey_flow_graph(
            b_edge_df, target_bait, target_prey, multi_flow_graph, mf_capacity, node_map
        )

        ssoe = multi_flow_graph.add_edge(super_source, source)
        ssie = multi_flow_graph.add_edge(sink, super_sink)
        mf_capacity[ssoe] = smallest_edge_count
        mf_capacity[ssie] = smallest_edge_count

        network_data[net_i] = {
            "source": source,
            "sink": sink,
            "flow_node_map_index": flow_node_map_index,
            "flow_node_map": flow_node_map,
            "node_map": node_map,
            "node_map_idx": node_map_idx,
        }

    return (
        multi_flow_graph,
        network_data,
        super_source,
        super_sink,
        mf_capacity,
        smallest_edge_count,
    )


def alternating_maxflow_multi_network(
    edge_list_a,
    edge_list_b,
    edge_deviation_threshold=0.1,
    min_flow=0.9,
):
    multi_flow_graph = [None, None]
    network_data = [None, None]
    super_source = [None, None]
    super_sink = [None, None]
    capacity = [None, None]
    smallest_edge_count = [None, None]
    selected_edges = [edge_list_a, [None] * len(edge_list_a)]

    (
        multi_flow_graph[0],
        network_data[0],
        super_source[0],
        super_sink[0],
        capacity[0],
        smallest_edge_count[0],
    ) = build_multi_network_flow_graph(edge_list_a, edge_list_b)

    i = 0
    current_least_flow = 0.0
    current_edge_deviation = 1.0

    while (
        current_least_flow < min_flow
        or current_edge_deviation > edge_deviation_threshold
    ):
        slot = i % 2
        next_slot = (i + 1) % 2

        residual = max_flow(
            multi_flow_graph[slot], super_source[slot], super_sink[slot], capacity[slot]
        )

        current_least_flow = 1
        current_edge_deviation = 1

        current_least_flow = min(
            [
                (capacity[slot][e] - residual[e]) / capacity[slot][e]
                for e in super_source[slot].out_edges()
            ]
        )

        single_network_degree = []
        per_subnetwork_edge_count = [
            capacity[slot][e] - residual[e] for e in super_source[slot].out_edges()
        ]
        for net_i in range(len(edge_list_a)):
            single_network_degree.append(
                sum(
                    capacity[slot][e]
                    for e in network_data[slot][net_i]["source"].out_edges()
                )
            )
        sample_balances = [
            abs(e / d - 1)
            for d, e in zip(single_network_degree, per_subnetwork_edge_count)
        ]
        current_edge_deviation = max(sample_balances)

        for net_i in range(len(edge_list_a)):
            selected_edges[next_slot][net_i] = graph_to_edge_df(
                extract_selected_edges(
                    multi_flow_graph[slot],
                    capacity[slot],
                    residual,
                    network_data[slot][net_i]["flow_node_map_index"],
                    network_data[slot][net_i]["node_map"],
                    network_data[slot][net_i]["source"],
                    network_data[slot][net_i]["sink"],
                ),
                network_data[slot][net_i]["node_map_idx"],
            )
            selected_edges[slot][net_i], selected_edges[next_slot][net_i] = (
                drop_exclusive_nodes(
                    selected_edges[slot][net_i], selected_edges[next_slot][net_i]
                )
            )

        (
            multi_flow_graph[next_slot],
            network_data[next_slot],
            super_source[next_slot],
            super_sink[next_slot],
            capacity[next_slot],
            smallest_edge_count[next_slot],
        ) = build_multi_network_flow_graph(
            selected_edges[next_slot], selected_edges[slot]
        )
        i += 1

    return [selected_edges[0][net_i] for net_i in range(len(edge_list_a))], [
        selected_edges[1][net_i] for net_i in range(len(edge_list_a))
    ]


def drop_exclusive_nodes(edge_df_a, edge_df_b):
    while True:
        if edge_df_a.empty or edge_df_b.empty:
            raise ValueError("One of the edge DataFrames is empty. All edges dropped.")
        n_edges = edge_df_a.shape[0] + edge_df_b.shape[0]
        all_baits = set(edge_df_a["bait"]) & set(edge_df_b["bait"])
        all_prey = set(edge_df_a["prey"]) & set(edge_df_b["prey"])

        edge_df_a = edge_df_a[
            (edge_df_a["bait"].isin(all_baits)) & (edge_df_a["prey"].isin(all_prey))
        ]
        edge_df_b = edge_df_b[
            (edge_df_b["bait"].isin(all_baits)) & (edge_df_b["prey"].isin(all_prey))
        ]
        if edge_df_a.shape[0] + edge_df_b.shape[0] == n_edges:
            break

    return edge_df_a, edge_df_b


def get_edge_list(
    full_detection_df,
    positive_limits,
    negative_limits,
    nodes_to_exclude,
):
    edge_list_pos = []
    edge_list_neg = []
    pair_order = []
    pos_edges = full_detection_df[
        ~full_detection_df["bait"].isin(nodes_to_exclude)
        & ~full_detection_df["prey"].isin(nodes_to_exclude)
    ]
    neg_edges = full_detection_df[
        ~full_detection_df["bait"].isin(nodes_to_exclude)
        & ~full_detection_df["prey"].isin(nodes_to_exclude)
    ]

    for pos_limit in positive_limits:
        for neg_limit in negative_limits:
            if pos_limit == "all":
                c_pos = pos_edges[pos_edges["n_observed"] != 0]
            else:
                c_pos = pos_edges[pos_edges["n_observed"] >= float(pos_limit)]

            c_neg = neg_edges[
                (neg_edges["n_tested"] >= int(neg_limit))
                & (neg_edges["n_observed"] == 0)
            ]

            pair_order.append((pos_limit, neg_limit))
            c_pos, c_neg = drop_exclusive_nodes(c_pos, c_neg)
            edge_list_pos.append(c_pos)
            edge_list_neg.append(c_neg)

    return edge_list_pos, edge_list_neg, pair_order


def main():
    test_set = snakemake.input.test_set
    validation_set = snakemake.input.validation_set

    full_detection_df = pd.read_parquet(snakemake.input.full_detection).rename(
        {"gene_name_bait": "bait", "gene_name_prey": "prey"}, axis=1
    )
    
    positive_limits = snakemake.params.positive_limits
    negative_limits = snakemake.params.negative_limits

    test_df = pd.read_csv(test_set, header=None, names=["bait", "prey"])
    validation_df = pd.read_csv(validation_set, header=None, names=["bait", "prey"])
    nodes_to_exclude = (
        set(test_df["bait"])
        | set(test_df["prey"])
        | set(validation_df["bait"])
        | set(validation_df["prey"])
    )

    edge_list_pos, edge_list_neg, pair = get_edge_list(
        full_detection_df,
        positive_limits,
        negative_limits,
        nodes_to_exclude,
    )

    selected_pos, selected_neg = alternating_maxflow_multi_network(
        edge_list_pos, edge_list_neg
    )

    balanced_edges_positive = snakemake.output.balanced_edges_positive
    balanced_edges_negative = snakemake.output.balanced_edges_negative
    for i, (output_pos, output_neg) in enumerate(
        zip(balanced_edges_positive, balanced_edges_negative)
    ):
        selected_pos[i].to_csv(output_pos, index=False, sep="\t")
        selected_neg[i].to_csv(output_neg, index=False, sep="\t")
