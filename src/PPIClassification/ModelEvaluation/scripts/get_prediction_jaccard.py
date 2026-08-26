import pandas as pd


def jaccard(set_a, set_b):
    union = set_a | set_b
    if not union:
        return float("nan")
    return len(set_a & set_b) / len(union)


def correct_pairs(df_subset, model_col):
    y_true = (df_subset["dataset"] == "Interaction").astype(int)
    y_pred = (df_subset[model_col] >= 0.5).astype(int)
    is_correct = y_pred == y_true
    pairs = zip(df_subset["bait"], df_subset["prey"])
    return {pair for pair, correct in zip(pairs, is_correct) if correct}


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
            "pos_limit": snakemake.wildcards.pos_limit,
            "neg_limit": snakemake.wildcards.neg_limit,
            "permutation": permutation,
            "label": label,
            "jaccard": jaccard(correct_pairs(sub_df, "HRNI"), correct_pairs(sub_df, "NO")),
            "n_pairs": len(sub_df),
        })

pd.DataFrame(rows).to_csv(snakemake.output.jaccard, sep="\t", index=False)
