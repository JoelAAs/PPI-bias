import pandas as pd
from sample_balance_multi_network_functions import degree_balace_edges


def main():
    pos_df = pd.read_csv(snakemake.input.balanced_positive_edges, sep="\t", dtype={"bait": "string", "prey": "string"})
    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")
    
    all_nodes = list(set(pos_df["bait"]) | set(pos_df["prey"]))
    pos_set = set(zip(pos_df["bait"], pos_df["prey"]))
    comp_df = pd.DataFrame(
        [(u, v) for u in all_nodes for v in all_nodes if u != v and (u, v) not in pos_set],
        columns=["bait", "prey"],
    )
    if not directed:
        comp_df["id"] = comp_df.apply(lambda row: "-".join(sorted(row)))
        comp_df.drop_duplicates(inplace=True)
        del comp_df["id"] # no duplicates for faster loading later

    _, selected_neg, _ = degree_balace_edges(pos_df, comp_df, min_flow=0.95, directed=directed)

    selected_neg.to_csv(snakemake.output.random_negative_edges, sep="\t", index=False)


if __name__ == "__main__":
    main()
