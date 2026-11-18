import sys
import numpy as np
import pandas as pd
from define_test_set import define_test_set
from sample_balance_multi_network_functions import sanity_check


def read_edges(path, directed):
    """
    Read a bait/prey edge list, drop self-pairs and duplicates. Undirected pairs are ordered bait < prey
    """
    df = pd.read_csv(path, sep="\t", dtype=str)[["bait", "prey"]]
    df = df[df["bait"] != df["prey"]]
    if not directed:
        swap = (df["bait"] > df["prey"]).to_numpy()
        df = pd.DataFrame({
            "bait": np.where(swap, df["prey"], df["bait"]),
            "prey": np.where(swap, df["bait"], df["prey"])
        })
    return df.drop_duplicates(ignore_index=True)


def drop_shared_pairs(pos_df, neg_df):
    """
    Pairs in both sets have conflicting labels, remove them from both
    """
    pos_idx = pd.MultiIndex.from_frame(pos_df)
    neg_idx = pd.MultiIndex.from_frame(neg_df)
    shared = pos_idx.intersection(neg_idx)
    print(f"Pairs in both positive and negative set (removed): {len(shared)}", flush=True)
    return pos_df[~pos_idx.isin(shared)], neg_df[~neg_idx.isin(shared)]


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    if snakemake.wildcards.network_type == "directional":
        directed=True
    else:
        directed=False

    pos_df = read_edges(snakemake.input.pos, directed)
    neg_df = read_edges(snakemake.input.neg, directed)
    
    pos_prots = set(pos_df["bait"]) | set(pos_df["prey"])
    neg_df = neg_df[(neg_df["bait"].isin(pos_prots)) & (neg_df["prey"].isin(pos_prots))]
    
    pos_df, neg_df = drop_shared_pairs(pos_df, neg_df)
    balanced_pos, balanced_neg = define_test_set(pos_df, neg_df, directed=directed) # is just the balancing, look over naming
    sanity_check([balanced_pos], [balanced_neg], [pos_df], [neg_df], directed=directed)
    print(f"Balanced positive: {balanced_pos.shape[0]}, negative: {balanced_neg.shape[0]}", flush=True)

    balanced = pd.concat(
        [balanced_pos.assign(interaction=True), balanced_neg.assign(interaction=False)],
        ignore_index=True
    )
    balanced[["bait", "prey", "interaction"]].to_csv(snakemake.output.balanced, sep="\t", index=False)


if __name__ == "__main__":
    main()
