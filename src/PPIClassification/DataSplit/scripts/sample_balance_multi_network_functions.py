import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def single_alternating_maxflow(positive_edge_df, negative_edge_df, min_flow, max_edges=None, seed=None):
    all_nodes = (
        set(positive_edge_df["bait"]) | set(positive_edge_df["prey"]) |
        set(negative_edge_df["bait"]) | set(negative_edge_df["prey"])
    )
    node_map, node_map_idx = get_node_map(all_nodes)
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

    edge_count_msg_count = positive_edge_df.shape[0]

    if seed is not None:
        np.random.seed(seed)

    while percent_flow_value < min_flow or (
        max_edges is not None and any(se.shape[0] > max_edges for se in subset_edges)
    ):
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
        if edge_count_msg_count - subset_edges[0].shape[0] >= 1000:
            print(f"Flow value {percent_flow_value}: positive edges {subset_edges[0].shape[0]}, negative edges {subset_edges[1].shape[0]}, max edges {max_edges}")
            edge_count_msg_count = subset_edges[0].shape[0]

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
        if (
            max_edges is not None
            and any(se.shape[0] > max_edges for se in subset_edges)
            and percent_flow_value >= 0.999
        ):
            cap = current_capacity[(i + 1) % 2]
            nonzero_indices = np.where(cap.a > 0)[0]
            if len(nonzero_indices) > 0:
                cap.a[np.random.choice(nonzero_indices)] = 0
        i += 1

    pos, neg = drop_exclusive_nodes(subset_edges[0], subset_edges[1])
    return pos, neg



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
    pos_edges = pos_edges[pos_edges["n_observed"] != 0]
    neg_edges = full_detection_df[
        ~full_detection_df["bait"].isin(nodes_to_exclude)
        & ~full_detection_df["prey"].isin(nodes_to_exclude)
    ]
    neg_edges = neg_edges[neg_edges["n_observed"] == 0]

    for pos_limit in positive_limits:
        if not pos_limit == "all":
            c_pos = pos_edges[pos_edges["lower_bound_pod"] >= float(pos_limit)]
        else:
            c_pos = pos_edges

        for neg_limit in negative_limits:
            c_neg = neg_edges[
                (neg_edges["n_tested"] >= int(neg_limit))
                & (neg_edges["n_observed"] == 0)
            ]

            pair_order.append((pos_limit, neg_limit))
            b_pos, b_neg = drop_exclusive_nodes(c_pos, c_neg)
            print(f"Positive edges {b_pos.shape[0]}, negative edges {b_neg.shape[0]}")
            edge_list_pos.append(b_pos)
            edge_list_neg.append(b_neg)

    return edge_list_pos, edge_list_neg, pair_order


def _run_single_maxflow(args):
    openmp_set_num_threads(1)
    pos_df, neg_df, min_flow, max_edges = args
    return single_alternating_maxflow(pos_df, neg_df, min_flow, max_edges)


def parallel_maxflow_multi_network(
    edge_list_a,
    edge_list_b,
    edge_deviation_threshold=0.05,
    min_flow=0.95,
    n_workers=None,
):
    i = 0
    current_edge_deviations = np.ones(len(edge_list_a))

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        selected_networks = list(executor.map(
            _run_single_maxflow,
            [(a, b, min_flow, None) for a, b in zip(edge_list_a, edge_list_b)],
        ))

        while any(current_edge_deviations > edge_deviation_threshold):
            current_edges = np.array([pos.shape[0] for pos, _ in selected_networks])

            if 0 in current_edges:
                raise ValueError("One of the selected networks has zero edges. Cannot proceed.")

            current_min_edges = int(min(current_edges))
            current_edge_deviations = np.array([
                (pos.shape[0] - current_min_edges) / current_min_edges
                for pos, _ in selected_networks
            ])

            c_max_edges = current_min_edges * (1 + edge_deviation_threshold)
            futures = {
                executor.submit(_run_single_maxflow, (
                    selected_networks[net_i][0],
                    selected_networks[net_i][1],
                    min_flow,
                    c_max_edges,
                )): net_i
                for net_i, dev in enumerate(current_edge_deviations)
                if dev > edge_deviation_threshold
            }
            for fut in as_completed(futures):
                selected_networks[futures[fut]] = fut.result()

            msg = f"Iteration {i}\tEdges\tRebalance\tDeviation\tMax Edges: {current_min_edges}\n"
            for net_i, (pos, neg) in enumerate(selected_networks):
                msg += f"Network {net_i}\t{pos.shape[0]}\t{current_edge_deviations[net_i] > edge_deviation_threshold}\t{current_edge_deviations[net_i]:.2f}\n"
            print(msg)
            i += 1

    return (
        [pos for pos, _ in selected_networks],
        [neg for _, neg in selected_networks],
    )
     
     
def sanity_check(selected_pos, selected_neg, edge_list_pos, edge_list_neg):
    def _make_sure_same_bait_prey(df1, df2):
        bait_same = set(df1["bait"]) == set(df2["bait"])
        prey_same = set(df1["prey"]) == set(df2["prey"])
        assert bait_same, "Bait sets do not match between selected and original edge lists."
        assert prey_same, "Prey sets do not match between selected and original edge lists."
        
        
    for i, (pos_df, neg_df, orig_pos_df, orig_neg_df) in enumerate(zip(selected_pos, selected_neg, edge_list_pos, edge_list_neg)):
        pos_edges_set = set(zip(pos_df["bait"], pos_df["prey"]))
        neg_edges_set = set(zip(neg_df["bait"], neg_df["prey"]))
        orig_pos_edges_set = set(zip(orig_pos_df["bait"], orig_pos_df["prey"]))
        orig_neg_edges_set = set(zip(orig_neg_df["bait"], orig_neg_df["prey"]))

        assert pos_edges_set.issubset(orig_pos_edges_set), f"Selected positive edges for network {i} are not a subset of original positive edges."
        assert neg_edges_set.issubset(orig_neg_edges_set), f"Selected negative edges for network {i} are not a subset of original negative edges."
        assert pos_edges_set.isdisjoint(neg_edges_set), f"Selected positive and negative edges for network {i} overlap."
        
        _make_sure_same_bait_prey(pos_df, neg_df)       

