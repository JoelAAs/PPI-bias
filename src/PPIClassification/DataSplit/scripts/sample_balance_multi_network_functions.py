import itertools
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from graph_tool.all import Graph, openmp_set_num_threads
from graph_tool.flow import push_relabel_max_flow as max_flow

def generate_graph(edge_df, node_map, directed):
    """Build a graph-tool Graph from a bait/prey edge DataFrame.

    Args:
        edge_df: DataFrame with columns "bait" and "prey".
        node_map: Dict mapping protein name -> integer vertex index.
        directed: If True, create a directed graph; otherwise undirected.

    Returns:
        graph_tool.Graph with vertices and edges added.
    """
    g = Graph(directed=directed)
    g.add_vertex(len(node_map))

    bait_src = edge_df["bait"].map(node_map).to_numpy()
    tar_prey = edge_df["prey"].map(node_map).to_numpy()
    edges = np.column_stack((bait_src, tar_prey)).astype(int)
    g.add_edge_list(edges)

    return g


def drop_exclusive_nodes(pos_edges, neg_edges, directed):
    """Remove proteins that appear in only one of the positive or negative edge sets.

    Iterates until convergence: at each step, retains only proteins present in
    both the positive and negative graphs, then filters edges accordingly.
    This ensures every protein has at least one positive and one negative edge,
    preventing degree-bias artefacts.

    Args:
        pos_edges: DataFrame of positive (HCI) bait/prey edges.
        neg_edges: DataFrame of negative (HCNI) bait/prey edges.

    Returns:
        Tuple (pos_edges, neg_edges) with exclusive-protein edges removed.
    """
    prev_edges = 0
    edges = pos_edges.shape[0] + neg_edges.shape[0]
    excluded_nodes = set()
    if not directed:
        original_proteins = (set(pos_edges["bait"]) | set(pos_edges["prey"])) & (set(neg_edges["bait"]) | set(neg_edges["prey"]))    
        while prev_edges != edges:
            prev_edges = edges
            common_proteins = (set(pos_edges["bait"]) | set(pos_edges["prey"])) & (set(neg_edges["bait"]) | set(neg_edges["prey"]))    
            pos_edges = pos_edges[(pos_edges["bait"].isin(common_proteins)) & (pos_edges["prey"].isin(common_proteins))]
            neg_edges = neg_edges[(neg_edges["bait"].isin(common_proteins)) & (neg_edges["prey"].isin(common_proteins))]
            edges = pos_edges.shape[0] + neg_edges.shape[0]
        excluded_nodes = original_proteins - common_proteins
        
    else:
        while prev_edges != edges:
            prev_edges = edges
            ex_bait = set(pos_edges["bait"]) ^ set(neg_edges["bait"])
            ex_prey = set(pos_edges["prey"]) ^ set(neg_edges["prey"])
            pos_edges = pos_edges[(~pos_edges["bait"].isin(ex_bait)) | (~pos_edges["prey"].isin(ex_prey))]
            neg_edges = neg_edges[(~neg_edges["bait"].isin(ex_bait)) | (~neg_edges["prey"].isin(ex_prey))]
            edges = pos_edges.shape[0] + neg_edges.shape[0]
            excluded_nodes |= ex_bait | ex_prey

    return pos_edges.copy(), neg_edges.copy(), excluded_nodes 


def sample_balance_multiple_networks(
    edge_list_pos,
    edge_list_neg,
    edge_deviation_threshold=0.05,
    min_flow=.95,
    n_workers=None,
    directed=True
):
    """Degree-balance multiple positive/negative network pairs so edge counts are comparable.

    First balances each network independently in parallel using `degree_balace_edges`.
    Then iteratively re-balances any network whose positive-edge count deviates from
    the cross-network minimum by more than `edge_deviation_threshold`, capping edges
    at the current minimum until all networks are within tolerance.

    Args:
        edge_list_pos: List of DataFrames, one positive edge set per network.
        edge_list_neg: List of DataFrames, one negative edge set per network.
        edge_deviation_threshold: Max allowed fractional deviation from the minimum
            edge count across networks (default 0.05 = 5%).
        min_flow: Minimum fraction of max-flow that must be achieved before
            accepting a balanced network (passed to `degree_balace_edges`).
        n_workers: Number of parallel worker processes (None = CPU count).
        directed: Whether to treat edges as directed bait->prey.

    Returns:
        Tuple (pos_list, neg_list) of balanced edge DataFrames in the same order
        as the inputs.
    """
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(degree_balace_edges, pos_edge, neg_edge, min_flow, directed)
            for pos_edge, neg_edge in zip(edge_list_pos, edge_list_neg)
        ]
        degree_balanced_networks = [f.result()[:2] for f in futures]

        current_min_edges = min(pos.shape[0] for pos, _ in degree_balanced_networks)
        current_edge_deviations = np.array([
            (pos.shape[0] - current_min_edges) / current_min_edges
            for pos, _ in degree_balanced_networks
        ])

        i = 0
        while any(current_edge_deviations > edge_deviation_threshold):
            run_mask = [c > edge_deviation_threshold for c in current_edge_deviations]
            reruned = {
                net_i: executor.submit(
                    degree_balace_edges,
                    degree_balanced_networks[net_i][0],
                    degree_balanced_networks[net_i][1],
                    min_flow,
                    directed,
                    current_min_edges,
                )
                for net_i, run in enumerate(run_mask) if run
            }

            for net_i, future in reruned.items():
                degree_balanced_networks[net_i] = future.result()[:2]

            current_min_edges = min(pos.shape[0] for pos, _ in degree_balanced_networks)
            current_edge_deviations = np.array([
                (pos.shape[0] - current_min_edges) / current_min_edges
                for pos, _ in degree_balanced_networks
            ])
            msg = f"Iteration {i}\tEdges\tRebalance\tDeviation\tMax Edges: {current_min_edges}\n"
            for net_i, dev in enumerate(current_edge_deviations):
                pos, _ = degree_balanced_networks[net_i]
                msg += f"Network {net_i}\t{pos.shape[0]}\t{dev > edge_deviation_threshold}\t{dev:.2f}\n"
            print(msg)
            i += 1
     
    return (
        [pos for pos, _ in degree_balanced_networks],
        [neg for _, neg in degree_balanced_networks],
    )
     
     
     
def get_node_map(all_nodes):
    """Create bidirectional mappings between protein names and integer vertex indices.

    Args:
        all_nodes: Iterable of protein identifiers (gene names or UniProt IDs).

    Returns:
        Tuple (node_map, node_map_idx) where node_map maps name->index and
        node_map_idx maps index->name.
    """
    all_nodes = list(all_nodes)
    return {gene: i for i, gene in enumerate(all_nodes)}, {
        i: gene for i, gene in enumerate(all_nodes)
    }


def graph_to_edge_df(g, node_map_idx):
    """Convert a graph-tool Graph back to a bait/prey edge DataFrame.

    Args:
        g: graph_tool.Graph whose edges should be exported.
        node_map_idx: Dict mapping integer vertex index -> protein name.

    Returns:
        DataFrame with columns ["bait", "prey"].
    """
    # order doesn't matter for undirected and columns stay as bait prey
    edge_list = []
    for e in g.edges():
        bait_idx = int(e.source())
        prey_idx = int(e.target())
        bait = node_map_idx[bait_idx]
        prey = node_map_idx[prey_idx]
        edge_list.append((bait, prey))
    return pd.DataFrame(edge_list, columns=["bait", "prey"])


def build_directed_flow_graph(
    edge_df, target_bait, target_prey, g, capacity, node_map
):
    """Extend a graph-tool Graph into a directed max-flow network for degree balancing.

    For each bait, adds a source->bait_copy edge with capacity equal to the target
    bait out-degree. For each prey, adds a prey_copy->sink edge with capacity equal
    to the target prey in-degree. Each original bait->prey edge becomes a capacity-1
    edge between the corresponding copies, so max-flow selects the largest subset of
    negative edges that respects the positive-graph degree sequence.

    Args:
        edge_df: Negative edge candidates as a bait/prey DataFrame.
        target_bait: Array of target out-degrees, one per bait in `node_map` order.
        target_prey: Array of target in-degrees, one per prey in `node_map` order.
        g: graph_tool.Graph to extend (modified in-place).
        capacity: Edge property map for capacities (modified in-place).
        node_map: Dict mapping protein name -> original vertex index.

    Returns:
        Tuple (g, capacity, source, sink, flow_node_map_index, flow_node_map).
    """
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


def build_undirected_flow_graph(
    edge_df, target_degree, g, capacity, node_map
):
    """Extend a graph-tool Graph into a directed max-flow network for undirected degree balancing.

    Each undirected edge A-B is represented as two directed edges (A->B and B->A) with
    capacity 1. Each node is split into two copies: a "source-side" copy with capacity
    from source, and a "sink-side" copy with capacity to sink, both capped at the target
    degree. Max-flow then selects the largest subset of negative edges that matches the
    undirected degree sequence of the positive graph.

    Args:
        edge_df: Negative edge candidates as a bait/prey DataFrame.
        target_degree: Array of target degrees, one per node in `node_map` order.
        g: graph_tool.Graph to extend (modified in-place).
        capacity: Edge property map for capacities (modified in-place).
        node_map: Dict mapping protein name -> original vertex index.

    Returns:
        Tuple (g, capacity, source, sink, flow_node_map_index, flow_node_map).
    """
    
    if edge_df.empty:
        raise ValueError("Edge DataFrame is empty. Cannot build flow graph.")
    if any(edge_df["bait"].isin(node_map) == False) or any(
        edge_df["prey"].isin(node_map) == False
    ):
        raise ValueError("Some baits or preys in edge_df are not in node_map.")
    num_vertices = g.num_vertices()

    # 0: A->B and 1: B->A for split edges
    flow_node_map = {(0, i): (i + num_vertices) for i, _ in enumerate(target_degree)}
    flow_node_map.update(
        {
            (1, i): (i + num_vertices + len(target_degree))
            for i, _ in enumerate(target_degree)
        }
    )
    flow_node_map_index = {value: key for key, value in flow_node_map.items()}

    g.add_vertex(len(flow_node_map))
    source = g.add_vertex()
    sink = g.add_vertex()

    for node, cap in enumerate(target_degree):
        source_v = flow_node_map[(0, node)]
        sink_v   = flow_node_map[(1, node)]
        e = g.add_edge(source, source_v)
        capacity[e] = int(cap)
        
        e = g.add_edge(sink_v, sink)
        capacity[e] = int(cap)


    num_existing = g.num_edges()
    edges = np.zeros(shape=(edge_df.shape[0]*2, 2), dtype=int)  # *2 for A->B and B->A
    from_nodes = edge_df["bait"].map(node_map).to_numpy(dtype=int)
    to_nodes   = edge_df["prey"].map(node_map).to_numpy(dtype=int)

    for i, (node_a, node_b) in enumerate(zip(from_nodes, to_nodes)):
        # (0, A) -> (1, B)
        edges[i*2]     = [flow_node_map[(0, node_a)], flow_node_map[(1, node_b)]]
        # (0, B) -> (1, A)
        edges[i*2 + 1] = [flow_node_map[(0, node_b)], flow_node_map[(1, node_a)]]
        
    
    g.add_edge_list(edges)

    # Set capacity=1 for all newly added edges in bulk
    capacity.a[num_existing : num_existing + edge_df.shape[0]*2] = 1

    return g, capacity, source, sink, flow_node_map_index, flow_node_map


def extract_selected_edges( 
    flow_g, capacity, residual, flow_node_map_index, node_map, directed
): # TODO: check if this needs to be modified for directed graphs (select greedily half edges are selected as edges)
    g = Graph(directed=directed)
    g.add_vertex(len(node_map))

    n_verts = flow_g.num_vertices()
    valid = np.zeros(n_verts, dtype=bool)
    orig_idx = np.zeros(n_verts, dtype=int)
    vids = np.array(list(flow_node_map_index.keys()))  # edge_id
    idxs = np.array([v[1] for v in flow_node_map_index.values()])  # (0, nid) or (1, nid)
    valid[vids] = True
    orig_idx[vids] = idxs

    edges = flow_g.get_edges(
        [flow_g.edge_index]
    )  # (n_edges, 3): [source, target, edge_idx]
    edge_flows = capacity.a[edges[:, 2]] - residual.a[edges[:, 2]]

    mask = (edge_flows == 1) & valid[edges[:, 0]] & valid[edges[:, 1]] # Will add all edges greedily for unidirected
    selected = orig_idx[edges[mask, :2]]

    if len(selected):
        g.add_edge_list(selected)
    return g

def get_degree(g):
    """Return the out-degree and in-degree arrays for all vertices in a directed graph.

    Args:
        g: Directed graph_tool.Graph.

    Returns:
        Tuple (out_degrees, in_degrees) as int64 numpy arrays indexed by vertex id.
    """
    return g.get_out_degrees(g.get_vertices()).astype(np.int64), g.get_in_degrees(
        g.get_vertices()
    ).astype(np.int64)
    

def degree_balace_edges(pos_edges, neg_edges, min_flow, directed, max_edges=None, seed=None):
    """Degree-balance a single positive/negative edge pair using alternating max-flow.

    Builds a flow network where the positive graph's per-node degree sequence acts as
    capacity constraints on the negative edges. Alternates between balancing negative
    edges against the positive degree sequence and positive edges against the negative
    degree sequence until `min_flow` fraction of max-flow is achieved. If `max_edges`
    is set, also randomly zeroes a capacity whenever the edge count exceeds the cap.

    Args:
        pos_edges: DataFrame of positive (HCI) bait/prey edges.
        neg_edges: DataFrame of negative (HCNI) bait/prey edges.
        min_flow: Minimum fraction of max-flow that must be reached before stopping
            (e.g. 0.95 means ≥95% of the theoretical maximum flow).
        directed: If True, balance bait out-degrees and prey in-degrees separately;
            if False, balance the combined undirected degree of each node.
        max_edges: Optional cap on the number of positive (or negative) edges in
            the final result. When set, the loop continues until both the flow
            threshold and the edge-count cap are satisfied.
        seed: Optional numpy random seed for reproducibility when randomly zeroing
            capacities.

    Returns:
        Tuple (pos_edges, neg_edges) of degree-balanced DataFrames with exclusive
        proteins removed (via `drop_exclusive_nodes`).
    """
    all_nodes = (
        set(pos_edges["bait"]) | set(pos_edges["prey"])
        | set(neg_edges["bait"]) | set(neg_edges["prey"])
    )
    node_map, node_map_idx = get_node_map(all_nodes)
    g_pos = generate_graph(pos_edges, node_map, directed)
    
    if directed:
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
        ) = build_directed_flow_graph(
            edge_df=neg_edges,
            target_bait=target_bait_pos,
            target_prey=target_prey_pos,
            g=mf_g_pos,
            capacity=capacity_pos,
            node_map=node_map,
        )
    else:
        target_degree = g_pos.get_total_degrees(g_pos.get_vertices()).astype(np.int64)
        mf_g_pos = Graph(directed=True)
        capacity_pos = mf_g_pos.new_edge_property("int")
        (
            mf_g_pos,
            capacity_pos,
            source_pos,
            sink_pos,
            flow_node_map_index_pos,
            flow_node_map_pos,
        ) = build_undirected_flow_graph(
            edge_df=neg_edges,
            target_degree=target_degree,
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
    subset_edges = [pos_edges, neg_edges]

    edge_count_msg_count = pos_edges.shape[0]

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
            directed=directed,
        )
        g_next = Graph(directed=directed)
        capacity_next = g_next.new_edge_property("int")
        subset_edges[(i + 1) % 2] = graph_to_edge_df(g_selected, node_map_idx)
        if edge_count_msg_count - subset_edges[0].shape[0] >= 1000:
            print(f"Flow value {percent_flow_value}: positive edges {subset_edges[0].shape[0]}, negative edges {subset_edges[1].shape[0]}, max edges {max_edges}")
            edge_count_msg_count = subset_edges[0].shape[0]
        if directed:
            new_target_bait, new_target_prey = get_degree(g_selected)
            (
                current_g[(i + 1) % 2],
                current_capacity[(i + 1) % 2],
                current_source[(i + 1) % 2],
                current_sink[(i + 1) % 2],
                current_flow_map_index[(i + 1) % 2],
                current_flow_map[(i + 1) % 2],
            ) = build_directed_flow_graph(
                edge_df=subset_edges[i % 2],
                target_bait=new_target_bait,
                target_prey=new_target_prey,
                g=g_next,
                capacity=capacity_next,
                node_map=node_map,
            )
        else:
            new_target_degree = g_selected.get_total_degrees(g_selected.get_vertices()).astype(np.int64)
            (
                current_g[(i + 1) % 2],
                current_capacity[(i + 1) % 2],
                current_source[(i + 1) % 2],
                current_sink[(i + 1) % 2],
                current_flow_map_index[(i + 1) % 2],
                current_flow_map[(i + 1) % 2],
            ) = build_undirected_flow_graph(
                edge_df=subset_edges[i % 2],
                target_degree=new_target_degree,
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

    pos, neg, excluded = drop_exclusive_nodes(subset_edges[0], subset_edges[1], directed)
    return pos, neg, excluded


def sanity_check(selected_pos, selected_neg, edge_list_pos, edge_list_neg, directed=True):
    """Assert structural invariants on balanced edge sets against their originals.

    Checks per network:
    - Selected positive edges are a subset of original positive edges.
    - Selected negative edges are a subset of original negative edges.
    - Selected positive and negative edges are disjoint.
    - The node sets (directed: bait and prey separately; undirected: all nodes combined)
      match between the selected positive and negative DataFrames.

    Args:
        selected_pos: List of balanced positive edge DataFrames.
        selected_neg: List of balanced negative edge DataFrames.
        edge_list_pos: List of original positive edge DataFrames (same order).
        edge_list_neg: List of original negative edge DataFrames (same order).
        directed: If True, treat (A, B) and (B, A) as distinct edges and check bait/prey
            sets separately. If False, treat edges as unordered and check the union node set.

    Raises:
        AssertionError: On the first violated invariant, with a message naming the network.
    """
    def _edge_set(df):
        if directed:
            return set(zip(df["bait"], df["prey"]))
        return {frozenset((a, b)) for a, b in zip(df["bait"], df["prey"])}

    def _make_sure_same_nodes(df1, df2, net_i):
        if directed:
            assert set(df1["bait"]) == set(df2["bait"]), \
                f"Bait sets do not match for network {net_i}."
            assert set(df1["prey"]) == set(df2["prey"]), \
                f"Prey sets do not match for network {net_i}."
        else:
            nodes1 = set(df1["bait"]) | set(df1["prey"])
            nodes2 = set(df2["bait"]) | set(df2["prey"])
            assert nodes1 == nodes2, \
                f"Node sets do not match for network {net_i}."

    for i, (pos_df, neg_df, orig_pos_df, orig_neg_df) in enumerate(
        zip(selected_pos, selected_neg, edge_list_pos, edge_list_neg)
    ):
        pos_edges_set      = _edge_set(pos_df)
        neg_edges_set      = _edge_set(neg_df)
        orig_pos_edges_set = _edge_set(orig_pos_df)
        orig_neg_edges_set = _edge_set(orig_neg_df)

        assert pos_edges_set.issubset(orig_pos_edges_set), \
            f"Selected positive edges for network {i} are not a subset of original positive edges."
        assert neg_edges_set.issubset(orig_neg_edges_set), \
            f"Selected negative edges for network {i} are not a subset of original negative edges."
        assert pos_edges_set.isdisjoint(neg_edges_set), \
            f"Selected positive and negative edges for network {i} overlap."

        _make_sure_same_nodes(pos_df, neg_df, i)


def get_edge_list(
    full_detection_df,
    positive_limits,
    negative_limits,
    nodes_to_exclude,
):
    """Build positive and negative edge lists for every combination of confidence thresholds.

    Iterates over the Cartesian product of `positive_limits` × `negative_limits`. For
    each pair the detection DataFrame is filtered into a positive and a negative set, with
    any protein in `nodes_to_exclude` removed from both (this keeps test/validation
    proteins out of training).

    Positive filter:
        - ``pos_limit == "all"``: any pair with ``n_observed >= 1`` (no POD threshold).
        - numeric ``pos_limit``: ``lower_bound_pod >= pos_limit``.

    Negative filter (always): ``n_tested >= neg_limit`` AND ``n_observed == 0``.

    Positives and negatives are disjoint by construction because positives require
    ``n_observed >= 1`` and negatives require ``n_observed == 0``.

    Args:
        full_detection_df: DataFrame with columns ``"bait"``, ``"prey"``,
            ``"lower_bound_pod"``, ``"n_tested"``, ``"n_observed"`` (gene-name or
            UniProt columns should be renamed to ``"bait"``/``"prey"`` before calling).
        positive_limits: List of lower_bound_pod thresholds or the sentinel ``"all"``.
        negative_limits: List of integer n_tested thresholds.
        nodes_to_exclude: Set of protein identifiers to exclude from both edge sets
            (typically the union of bait and prey nodes in the test and validation sets).

    Returns:
        Tuple ``(edge_list_pos, edge_list_neg, pairs)`` where:
        - ``edge_list_pos``: list of positive edge DataFrames (columns: bait, prey).
        - ``edge_list_neg``: list of negative edge DataFrames (columns: bait, prey).
        - ``pairs``: list of ``(pos_limit, neg_limit)`` tuples in the same order,
          matching the Snakemake ``expand()`` output order.
    """
    in_excluded = (
        full_detection_df["bait"].isin(nodes_to_exclude)
        | full_detection_df["prey"].isin(nodes_to_exclude)
    )
    train_df = full_detection_df[~in_excluded]

    edge_list_pos = []
    edge_list_neg = []
    pairs = []

    for pos_limit, neg_limit in itertools.product(positive_limits, negative_limits):
        if pos_limit == "all":
            pos_mask = train_df["n_observed"] >= 1
        else:
            pos_mask = train_df["lower_bound_pod"] >= pos_limit

        neg_mask = (train_df["n_tested"] >= neg_limit) & (train_df["n_observed"] == 0)

        edge_list_pos.append(train_df.loc[pos_mask, ["bait", "prey"]].copy())
        edge_list_neg.append(train_df.loc[neg_mask, ["bait", "prey"]].copy())
        pairs.append((pos_limit, neg_limit))

    return edge_list_pos, edge_list_neg, pairs
