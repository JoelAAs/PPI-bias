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



def build_bait_prey_flow_graph(negative_df, target_bait, target_prey, g, node_map):
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
    #g = Graph(directed=True)
    num_edges = g.num_edges()
    
    flow_node_map = {(0,i):(i+ num_edges) for i, _ in enumerate(target_bait)}
    flow_node_map.update({(1,i):(i + num_edges+ len(target_bait)) for i, _ in enumerate(target_prey)})
    flow_node_map_index = {value:key for key, value in flow_node_map.items()}
    
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

    src_baits = negative_df["bait"].map(node_map).to_numpy()
    tar_prey = negative_df["prey"].map(node_map).to_numpy()

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


def extract_selected_edges(flow_g, capacity, residual, flow_node_map_index, node_map_index, source, sink):
    selected_edges = []

    for e in flow_g.edges():
        u = e.source()
        v = e.target()

        if u == source or v == sink:
            continue

        flow = capacity[e] - residual[e] # if not in residual then edge is chosen

        if flow == 1:
            bait = node_map_index[flow_node_map_index[int(u)][1]]
            prey = node_map_index[flow_node_map_index[int(v)][1]]
            selected_edges.append((bait, prey))

    return pd.DataFrame(selected_edges, columns=["bait", "prey"])

