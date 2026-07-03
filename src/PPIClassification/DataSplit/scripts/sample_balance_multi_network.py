import pandas as pd
from sample_balance_multi_network_functions import *

def main():
    # still called bait and prey in columns but they are protein A and protein B
    workers = snakemake.threads
    test_set = snakemake.input.test_set
    validation_set = snakemake.input.validation_set
    
    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError("{network_type} is not a valid network type.")

    full_detection_df = pd.read_parquet(snakemake.input.full_detection)
    bait_col = next(c for c in full_detection_df.columns if c.endswith("_bait"))
    prey_col = next(c for c in full_detection_df.columns if c.endswith("_prey"))
    full_detection_df = full_detection_df.rename({bait_col: "bait", prey_col: "prey"}, axis=1)

    positive_limits = snakemake.params.positive_limits
    negative_limits = snakemake.params.negative_limits

    test_df = pd.read_csv(test_set, header=None, names=["bait", "prey"], sep="\t")
    validation_df = pd.read_csv(
        validation_set, header=None, names=["bait", "prey"], sep="\t"
    )
    nodes_to_exclude = (
        set(test_df["bait"])
        | set(test_df["prey"])
        | set(validation_df["bait"])
        | set(validation_df["prey"])
    )

    edge_list_pos, edge_list_neg, pair = get_edge_list(
        full_detection_df,
        positive_limits,
        negative_limits,
        nodes_to_exclude,
    )

    log_file = snakemake.log[0]
    selected_pos, selected_neg = sample_balance_multiple_networks(
        edge_list_pos, edge_list_neg, log_file, n_workers=workers, directed=directed
    )

    sanity_check(selected_pos, selected_neg, edge_list_pos, edge_list_neg, directed=directed, log_file=log_file)
    
    balanced_edges_positive = snakemake.output.balanced_edges_positive
    balanced_edges_negative = snakemake.output.balanced_edges_negative
    for i, (output_pos, output_neg) in enumerate(
        zip(balanced_edges_positive, balanced_edges_negative)
    ):
        selected_pos[i].to_csv(output_pos, index=False, sep="\t")
        selected_neg[i].to_csv(output_neg, index=False, sep="\t")


if __name__ == "__main__":
    main()
