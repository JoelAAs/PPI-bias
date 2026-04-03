import pandas as pd
import argparse
import numpy as np
from scipy.stats import spearmanr
from graph_tool.all import Graph, openmp_set_num_threads
from graph_tool.flow import push_relabel_max_flow as max_flow

# target_degree is a list with degree where each position correlates to the node idx in the network
#
#
#
#

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
    return g.get_out_degrees(g.get_vertices()).astype(np.int64), g.get_in_degrees(g.get_vertices()).astype(np.int64)

def get_node_map(all_nodes):
    all_nodes = list(all_nodes)
    return {gene:i for i, gene in enumerate(all_nodes)}, {i:gene for i, gene in enumerate(all_nodes)}


def build_bait_prey_flow_graph(edge_df, target_bait, target_prey, g, node_map):
    """_summary_

    Args:
        negative_df (negative_edges): dataframe of non-interaction edges
        target_bait (list): degree per bait node
        target_prey (list): degree per prey node
        g (Graph): directed graph-tool graph

    Returns:
        _type_: _description_
    """
    # 0 is bait, 1 is prey in tuples
    # g = Graph(directed=True)
    num_edges = g.num_edges()

    flow_node_map = {(0, i): (i + num_edges) for i, _ in enumerate(target_bait)}
    flow_node_map.update(
        {(1, i): (i + num_edges + len(target_bait)) for i, _ in enumerate(target_prey)}
    )
    flow_node_map_index = {value: key for key, value in flow_node_map.items()}

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
    # Need to make sure that it has the same
    g = Graph(directed=True)
    g.add_vertex(len(node_map))

    selected_edges = []

    for e in flow_g.edges():
        u = e.source()
        v = e.target()

        if u == source or v == sink:
            continue

        flow = capacity[e] - residual[e]  # if not in residual then edge is chosen

        if flow == 1:
            bait = flow_node_map_index[int(u)][1]
            prey = flow_node_map_index[int(v)][1]
            selected_edges.append((bait, prey))
    g.add_edge_list(selected_edges)
    return g


def single_alternating_maxflow(positive_edge_df, negative_edge_df, min_flow):
    shared_baits = set(negative_edge_df["bait"]) & set(positive_edge_df["bait"])
    shared_prey = set(negative_edge_df["prey"]) & set(positive_edge_df["prey"])
    
    node_map, node_map_idx = get_node_map(shared_baits | shared_prey)
    g_pos = generate_graph(positive_edge_df, node_map)
    target_bait_pos, target_prey_pos = get_degree(g_pos)
    
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
        Graph(directed=True),
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
            current_g[i % 2], current_source[i % 2], current_sink[i % 2], current_capacity[i % 2])
        percent_flow_value = sum(
                current_capacity[i % 2][e] - residual[e]
                for e in current_source[i % 2].out_edges()
            )/sum(current_capacity[i % 2][e] for e in current_source[i % 2].out_edges())
        
        g_selected = extract_selected_edges(
            current_g[i % 2],
            current_capacity[i % 2],
            residual,
            current_flow_map_index[i % 2],
            node_map,
            current_source[i % 2],
            current_sink[i % 2])
        subset_networks[(i+1) %2] = g_selected.copy()
        new_target_bait, new_target_prey = get_degree(g_selected)
        if current_g[(i+1) % 2] is None:
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
                Graph(directed=True),
                node_map,
            )
        else:
            current_capacity[(i+1) % 2] = update_capacity(
                current_g[(i+1) % 2],
                current_capacity[(i+1) % 2],
                new_target_bait,
                new_target_prey,
                current_flow_map[(i+1) % 2],
                current_sink[(i+1) % 2],
                current_source[(i+1) % 2]
            )

    return subset_networks[0], subset_networks[1], node_map_idx


def build_multi_network_flow_graph(edge_list_a, edge_list_b):
    if len(edge_list_a) != len(edge_list_b):
        raise ValueError("Positive and negative edge lists must have the same length.")
    min_flow = 0.9
    
    smallest_edge_count = min([
        single_alternating_maxflow(c_positive_edge_df, c_negative_edge_df, min_flow)[0].num_edges()
        for c_positive_edge_df, c_negative_edge_df in zip(positive_edge_list, negative_edge_list)
    ])
    
    multi_flow_graph = Graph(directed=True)
    for i, c_positive_edge_df, c_negative_edge_df in enumerate(zip(positive_edge_list, negative_edge_list)):
        
        
    
    