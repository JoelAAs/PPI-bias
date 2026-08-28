import pandas as pd


def jaccard(set_a, set_b):
    union = set_a | set_b
    if not union:
        return float("nan")
    return len(set_a & set_b) / len(union)


def select_pairs(df_subset, model_col, pair_set):
    y_pred = (df_subset[model_col] >= 0.5).astype(int)
    pairs = zip(df_subset["bait"], df_subset["prey"])
    if pair_set == "correct":
        y_true = (df_subset["dataset"] == "Interaction").astype(int)
        keep = y_pred == y_true
    elif pair_set == "all":
        # "all" pairs the model calls an interaction, regardless of ground truth.
        keep = y_pred == 1
    else:
        raise ValueError(f"{pair_set} is not a valid pair_set.")
    return {pair for pair, k in zip(pairs, keep) if k}


pair_set = snakemake.wildcards.pair_set

rows = []
for pred_file in snakemake.input.predictions:
    permutation = pred_file.split("/permuted/")[1].split("/")[0]
    df = pd.read_csv(pred_file, sep="\t", dtype={"bait": "string", "prey": "string"})

    groups = {"global": df}
    groups.update(dict(tuple(df.groupby("dataset"))))

    for label, sub_df in groups.items():
        rows.append({
            "dataset": snakemake.wildcards.dataset,
            "network_type": snakemake.wildcards.network_type,
            "pair_set": pair_set,
            "pos_limit": snakemake.wildcards.pos_limit,
            "neg_limit": snakemake.wildcards.neg_limit,
            "permutation": permutation,
            "label": label,
            "jaccard": jaccard(select_pairs(sub_df, "HRNI", pair_set), select_pairs(sub_df, "NO", pair_set)),
            "n_pairs": len(sub_df),
        })

pd.DataFrame(rows).to_csv(snakemake.output.jaccard, sep="\t", index=False)
