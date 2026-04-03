import pandas as pd
import argparse
import numpy as np
from scipy.stats import spearmanr
from graph_tool.all import Graph, openmp_set_num_threads
from graph_tool.flow import push_relabel_max_flow as max_flow


def generate_graph(edge_df, node_map):
    g = Graph(directed=True)
    g.add_vertex(len(node_map))

    bait_src = edge_df["bait"].map(node_map).to_numpy()
    tar_prey = edge_df["prey"].map(node_map).to_numpy()
    edges = np.column_stack((bait_src, tar_prey))
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
    num_edges = g.num_edges()

    flow_node_map = {(0, i): (i + num_edges) for i, _ in enumerate(target_bait)}
    flow_node_map.update(
        {(1, i): (i + num_edges + len(target_bait)) for i, _ in enumerate(target_prey)}
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

    src_baits = edge_df["bait"].map(node_map).to_numpy()
    tar_prey = edge_df["prey"].map(node_map).to_numpy()

    bait_ids = [flow_node_map[(0, i)] for i in src_baits]
    prey_ids = [flow_node_map[(1, i)] for i in tar_prey]

    edges = np.column_stack((bait_ids, prey_ids))

    for u, v in edges:
        e = g.add_edge(u, v)
        capacity[e] = 1
    return g, capacity, source, sink, flow_node_map_index, flow_node_map


def update_capacity(g, capacity, target_bait, target_prey, flow_node_map, sink, source):
    for bait, cap in enumerate(target_bait):
        bait_v = flow_node_map[(0, bait)]
        e = g.get_edge(source, bait_v)
        capacity[e] = int(cap)

    for prey, cap in enumerate(target_prey):
        prey_v = flow_node_map[(1, prey)]
        e = g.get_edge(prey_v, sink)
        capacity[e] = int(cap)

    return capacity


def extract_selected_edges(
    flow_g, capacity, residual, flow_node_map_index, node_map, source, sink
):
    g = Graph(directed=True)
    g.add_vertex(len(node_map))

    selected_edges = []

    for e in flow_g.edges():
        u = int(e.source())
        v = int(e.target())

        if u not in flow_node_map_index or v not in flow_node_map_index:
            continue

        flow = capacity[e] - residual[e]

        if flow == 1:
            bait = flow_node_map_index[u][1]
            prey = flow_node_map_index[v][1]
            selected_edges.append((bait, prey))

    g.add_edge_list(selected_edges)
    return g



def single_alternating_maxflow(positive_edge_df, negative_edge_df, min_flow):
    shared_baits = set(negative_edge_df["bait"]) & set(positive_edge_df["bait"])
    shared_prey = set(negative_edge_df["prey"]) & set(positive_edge_df["prey"])

    node_map, node_map_idx = get_node_map(shared_baits | shared_prey)
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
        negative_edge_df,
        target_bait_pos,
        target_prey_pos,
        mf_g_pos,
        capacity_pos,
        node_map,
    )

    i = 0
    percent_flow_value = 0
    current_g = [mf_g_pos, None]
    current_capacity = [capacity_pos, None]
    current_source = [source_pos, None]
    current_sink = [sink_pos, None]
    current_flow_map_index = [flow_node_map_index_pos, None]
    current_flow_map = [flow_node_map_pos, None]
    subset_networks = [mf_g_pos.copy(), None]

    while percent_flow_value < min_flow:

        residual = max_flow(
            current_g[i % 2],
            current_source[i % 2],
            current_sink[i % 2],
            current_capacity[i % 2],
        )
        percent_flow_value = sum(
            current_capacity[i % 2][e] - residual[e]
            for e in current_source[i % 2].out_edges()
        ) / sum(current_capacity[i % 2][e] for e in current_source[i % 2].out_edges())

        g_selected = extract_selected_edges(
            current_g[i % 2],
            current_capacity[i % 2],
            residual,
            current_flow_map_index[i % 2],
            node_map,
            current_source[i % 2],
            current_sink[i % 2],
        )
        subset_networks[(i + 1) % 2] = g_selected.copy()
        new_target_bait, new_target_prey = get_degree(g_selected)
        if current_g[(i + 1) % 2] is None:
            g_neg = Graph(directed=True)
            capacity_neg = g_neg.new_edge_property("int")
            (
                current_g[1],
                current_capacity[1],
                current_source[1],
                current_sink[1],
                current_flow_map_index[1],
                current_flow_map[1],
            ) = build_bait_prey_flow_graph(
                positive_edge_df,
                new_target_bait,
                new_target_prey,
                g_neg,
                capacity_neg,
                node_map,
            )
        else:
            current_capacity[(i + 1) % 2] = update_capacity(
                current_g[(i + 1) % 2],
                current_capacity[(i + 1) % 2],
                new_target_bait,
                new_target_prey,
                current_flow_map[(i + 1) % 2],
                current_sink[(i + 1) % 2],
                current_source[(i + 1) % 2],
            )
        i += 1

    return subset_networks[0], subset_networks[1], node_map_idx


def build_multi_network_flow_graph(edge_list_a, edge_list_b):
    if len(edge_list_a) != len(edge_list_b):
        raise ValueError("Positive and negative edge lists must have the same length.")
    min_flow = 0.9

    smallest_edge_count = min(
        single_alternating_maxflow(a_edge_df, b_edge_df, min_flow)[0].num_edges()
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
            flow_node_map,   # captured here
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
    g_selected = [[None]*len(edge_list_a), [None]*len(edge_list_a)]
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

        for net_i in range(len(edge_list_a)):
            src = network_data[slot][net_i]["source"]

            percent_flow_value = sum(
                [capacity[slot][e] - residual[e] for e in src.out_edges()]
            ) / sum([capacity[slot][e] for e in src.out_edges()])
            if current_least_flow > percent_flow_value:
                current_least_flow = percent_flow_value

        per_subnetwork_edge_count = [
            capacity[slot][e] - residual[e] for e in super_source.out_edges()
        ]
        try:
            current_edge_deviation = (
                max(per_subnetwork_edge_count) / min(per_subnetwork_edge_count) - 1
            )
        except ZeroDivisionError:
            raise ValueError(
                "One of the sub-networks has zero selected edges, cannot compute edge deviation."
            )

        for net_i in range(len(edge_list_a)):
            g_selected[slot][i] = extract_selected_edges(
                multi_flow_graph[slot],
                capacity[slot],
                residual,
                network_data[slot][net_i]["flow_node_map_index"],
                network_data[slot][net_i]["node_map"],
                network_data[slot][net_i]["source"],
                network_data[slot][net_i]["sink"],
        )

        if multi_flow_graph[next_slot] is None:
            (
                multi_flow_graph[1],
                network_data[1],
                super_source[1],
                super_sink[1],
                capacity[1],
                smallest_edge_count[1],
            ) = build_multi_network_flow_graph(g_selected[slot][i], edge_list_a) # expects pandas df, not graph-tool graph; need to adjust build_multi_network_flow_graph to handle this case
        else:
            for net_i in range(n):
                capacity[next_slot] = update_capacity(
                    multi_flow_graph[next_slot],
                    capacity[next_slot],
                    *get_degree(g_selected[slot][net_i]),
                    network_data[next_slot][net_i]["flow_node_map"],
                    network_data[next_slot][net_i]["sink"],
                    network_data[next_slot][net_i]["source"],
                )

        i += 1

    return g_selected[0], g_selected[1]
