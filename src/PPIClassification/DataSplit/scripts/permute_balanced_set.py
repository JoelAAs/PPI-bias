import sys
import pandas as pd
from sample_balance_multi_network_functions import degree_balace_edges


def main():
    sys.stdout = open(snakemake.log[0], "w")
    sys.stderr = sys.stdout

    permutation = int(snakemake.wildcards.permutation)
    fraction = snakemake.params.fraction
    base_seed = snakemake.params.base_seed
    min_flow = snakemake.params.min_flow
    seed = base_seed + permutation

    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise CommandError("{network_type} is not a valid network type.")


    pos_df = pd.read_csv(snakemake.input.balanced_pos, sep="\t")
    neg_df = pd.read_csv(snakemake.input.balanced_neg, sep="\t")

    current_edges = pos_df.shape[0]
    target_edges = int(current_edges * fraction)
    pos_permut_df, neg_permut_df, discarded_nodes = degree_balace_edges(pos_df, neg_df, min_flow, max_edges=target_edges, directed=directed, seed=seed)

    pos_permut_df.to_csv(snakemake.output.permuted_pos, index=False, sep="\t")
    neg_permut_df.to_csv(snakemake.output.permuted_neg, index=False, sep="\t")


if __name__ == "__main__":
    main()
