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

edges = {
    permutation_index(path): edge_id(read_positive_edges(path), directed)
    for path in snakemake.input.edge_files
}

rows = []
for perm_a, perm_b in itertools.combinations(sorted(edges, key=int), 2):
    rows.append({
        "permutation_a": perm_a,
        "permutation_b": perm_b,
        "jaccard_index": jaccard_index(edges[perm_a], edges[perm_b]),
    })

pd.DataFrame(rows).to_csv(snakemake.output.jaccard_distances, sep="\t", index=False)
