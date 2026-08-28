import itertools
import pandas as pd


def permutation_index(path):
    return path.split("/permuted/")[1].split("/")[0]


def read_positive_edges(path):
    return pd.read_csv(path, sep="\t", dtype={"bait": "string", "prey": "string"})


def edge_id(df, directed):
    if directed:
        return df["bait"] + "-" + df["prey"]
    return df.apply(lambda row: "-".join(sorted([row["bait"], row["prey"]])), axis=1)


def jaccard_index(ids_a, ids_b):
    set_a, set_b = set(ids_a), set(ids_b)
    union = set_a | set_b
    return (len(set_a & set_b) / len(union)) if union else float("nan")


directed = snakemake.wildcards.network_type == "directional"
pos_limits = snakemake.params.pos_limits

edges_by_pos_limit = {}
for pos_limit in pos_limits:
    paths = [p for p in snakemake.input.permuted_pos if f"_poslim_{pos_limit}_pos.csv" in p]
    edges_by_pos_limit[str(pos_limit)] = {
        permutation_index(path): edge_id(read_positive_edges(path), directed)
        for path in paths
    }

rows = []
for pos_from, pos_to in zip(pos_limits, pos_limits[1:]):
    edges_from = edges_by_pos_limit[str(pos_from)]
    edges_to = edges_by_pos_limit[str(pos_to)]
    label = f"{pos_from}-{pos_to}"
    for perm_a, perm_b in itertools.product(sorted(edges_from, key=int), sorted(edges_to, key=int)):
        rows.append({
            "permutation_a": perm_a,
            "permutation_b": perm_b,
            "pos_limit": label,
            "jaccard_index": jaccard_index(edges_from[perm_a], edges_to[perm_b]),
        })

pd.DataFrame(rows).to_csv(snakemake.output.jaccard_distances, sep="\t", index=False)
