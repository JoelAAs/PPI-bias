import itertools
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr


def read_edges(path):
    return pd.read_csv(path, sep="\t", dtype={"bait": "string", "prey": "string"})


def degree_dicts(df, directed):
    graph_type = nx.DiGraph if directed else nx.Graph
    G = nx.from_pandas_edgelist(df, "bait", "prey", create_using=graph_type)
    if directed:
        return dict(G.out_degree()), dict(G.in_degree())
    # Undirected: "bait"/"prey" are just the two columns of an unordered edge, so a
    # protein's real degree counts both regardless of which column it landed in.
    deg = dict(G.degree())
    return deg, deg


def edge_set(df, directed):
    # Undirected: canonicalize (bait, prey) vs (prey, bait) to the same id, matching the
    # "A->B and B->A is considered the same" semantics - bait/prey column order isn't
    # otherwise guaranteed consistent across permutation files.
    pairs = zip(df["bait"], df["prey"])
    return set(pairs) if directed else set(tuple(sorted(p)) for p in pairs)


def jaccard_edges(df_a, df_b, directed):
    set_a, set_b = edge_set(df_a, directed), edge_set(df_b, directed)
    union = set_a | set_b
    return (len(set_a & set_b) / len(union)) if union else float("nan")


def spearman_degree(deg_a, deg_b):
    nodes = sorted(set(deg_a) | set(deg_b))
    if len(nodes) < 2:
        return float("nan")
    a = [deg_a.get(n, 0) for n in nodes]
    b = [deg_b.get(n, 0) for n in nodes]
    return spearmanr(a, b).correlation


def permutation_index(path):
    return path.split("/permuted/")[1].split("/")[0]


network_type = snakemake.wildcards.network_type
if network_type == "directional":
    directed = True
elif network_type == "undirectional":
    directed = False
else:
    raise ValueError(f"{network_type} is not a valid network type.")

groups = {
    "pos": snakemake.input.perm_pos,
    "neg_hrni": snakemake.input.perm_neg,
    "neg_random": snakemake.input.perm_random,
}

rows = []
for label, paths in groups.items():
    edges = {permutation_index(p): read_edges(p) for p in paths}
    degrees = {perm: degree_dicts(df, directed) for perm, df in edges.items()}
    for perm_a, perm_b in itertools.combinations(sorted(edges, key=int), 2):
        df_a, df_b = edges[perm_a], edges[perm_b]
        bait_a, prey_a = degrees[perm_a]
        bait_b, prey_b = degrees[perm_b]
        rows.append({
            "dataset": snakemake.wildcards.dataset,
            "network_type": snakemake.wildcards.network_type,
            "pos_limit": snakemake.wildcards.pos_limit,
            "neg_limit": snakemake.wildcards.neg_limit,
            "label": label,
            "permutation_a": perm_a,
            "permutation_b": perm_b,
            "jaccard_edges": jaccard_edges(df_a, df_b, directed),
            "spearman_bait": spearman_degree(bait_a, bait_b),
            "spearman_prey": spearman_degree(prey_a, prey_b),
            "n_a": len(df_a),
            "n_b": len(df_b),
        })

pd.DataFrame(rows).to_csv(snakemake.output.similarity, sep="\t", index=False)
