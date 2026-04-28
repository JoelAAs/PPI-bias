import pandas as pd
from sample_balance_multi_network import single_alternating_maxflow


def main():
    permutation = int(snakemake.wildcards.permutation)
    fraction = snakemake.params.fraction
    base_seed = snakemake.params.base_seed
    min_flow = snakemake.params.min_flow
    seed = base_seed + permutation

    pos_df = pd.read_csv(snakemake.input.balanced_pos, sep="\t")
    neg_df = pd.read_csv(snakemake.input.balanced_neg, sep="\t")

    current_edges = pos_df.shape[0]
    target_edges = int(current_edges * fraction)
    pos_permut_df, neg_permut_df = single_alternating_maxflow(pos_df, neg_df, min_flow, target_edges, seed=seed)

    pos_permut_df.to_csv(snakemake.output.permuted_pos, index=False, sep="\t")
    neg_permut_df.to_csv(snakemake.output.permuted_neg, index=False, sep="\t")

main()
