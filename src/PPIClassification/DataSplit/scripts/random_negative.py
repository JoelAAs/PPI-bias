import pandas as pd
import itertools
from sample_balance_multi_network_functions import degree_balace_edges


def get_no_edges(pos_df, directed, edges_to_avioid=None):
    all_nodes = list(set(pos_df["bait"]) | set(pos_df["prey"]))
    pos_set = set(zip(pos_df["bait"], pos_df["prey"])) 
    if edges_to_avioid:
        pos_set |= edges_to_avioid
    
    if not directed:
        pos_set |= {(b,a) for (a,b) in pos_set}
    
    pairs = itertools.permutations(all_nodes, 2) if directed else itertools.combinations(all_nodes, 2)
    comp_df = pd.DataFrame([p for p in pairs if p not in pos_set], columns=["bait", "prey"])

    _, selected_neg, _ = degree_balace_edges(pos_df, comp_df, min_flow=0.95, directed=directed)
    return selected_neg


def main():
    pos_df = pd.read_csv(snakemake.input.balanced_positive_edges, sep="\t", dtype={"bait": "string", "prey": "string"})
    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")
    
    selected_neg = get_no_edges(pos_df, directed)

    selected_neg.to_csv(snakemake.output.random_negative_edges, sep="\t", index=False)


if __name__ == "__main__":
    main()
