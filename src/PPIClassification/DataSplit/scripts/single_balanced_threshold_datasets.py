import pandas as pd
from sample_balance_multi_network_functions import *


def main():
    # still called bait and prey in columns but they can be protein A and protein B
    test_set = snakemake.input.test_set

    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError("{network_type} is not a valid network type.")

    full_detection_df = pd.read_parquet(snakemake.input.full_detection)

    source = snakemake.params.source
    pos_limit = snakemake.wildcards.pos_limit
    neg_limit = snakemake.wildcards.neg_limit

    test_df = pd.read_csv(test_set, header=None, names=["bait", "prey"], sep="\t")
    nodes_to_exclude = (
        set(test_df["bait"])
        | set(test_df["prey"])
    )

    edge_pos, edge_neg = get_edge_list(
        full_detection_df,
        source["positive"],
        source["negative"],
        pos_limit,
        neg_limit,
        nodes_to_exclude,
    )
    network_id = f"poslim_{pos_limit}_neglim_{neg_limit}"

    log_file = snakemake.log[0]
    selected_pos, selected_neg, _ = degree_balace_edges(
        edge_pos, edge_neg, 0.95, directed, None, None, log_file, network_id
    )

    sanity_check([selected_pos], [selected_neg], [edge_pos], [edge_neg], directed=directed, log_file=log_file)

    selected_pos.to_csv(snakemake.output.balanced_edges_positive, index=False, sep="\t")
    selected_neg.to_csv(snakemake.output.balanced_edges_negative, index=False, sep="\t")


if __name__ == "__main__":
    main()
