import re

import pandas as pd

MODEL_COLUMN_RE = re.compile(r"^train_(?P<train_type>hrni|no)_poslim(?P<pos_limit>all|[0-9.]+)_neglim(?P<neg_limit>[0-9.]+)$")

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.predictions, sep="\t", dtype={"bait": "string", "prey": "string"})
    model_columns = [c for c in df.columns if MODEL_COLUMN_RE.match(c)]

    train_type_columns = {}
    for column in model_columns:
        train_type = MODEL_COLUMN_RE.match(column)["train_type"]
        train_type_columns.setdefault(train_type, []).append(column)

    # "dataset" here is the ground-truth test-set label (Interaction/Negative_HRNI/Negative_NO)
    # written by evaluate_all_model_predictions.py, not the ms/y2h/flat dataset wildcard.
    is_interaction = df["dataset"] == "Interaction"
    is_noninteraction = df["dataset"].isin(["Negative_HRNI", "Negative_NO"])
    truth = is_interaction.astype(int)

    frames = []
    for train_type, columns in train_type_columns.items():
        n_models = len(columns)
        predicted_positive = (df[columns] >= 0.5).astype(int)
        # correct count per test pair, summed over this permutation's pos_limit x neg_limit
        # models (0..n_models) for the given train_type.
        n_correct_per_pair = predicted_positive.eq(truth, axis=0).sum(axis=1)

        inter_stats = (
            pd.DataFrame({"protein": df.loc[is_interaction, "bait"], "n_correct": n_correct_per_pair[is_interaction]})
            .groupby("protein")["n_correct"]
            .agg(n_correct_interaction="sum", n_pairs_interaction="count")
        )
        non_stats = (
            pd.DataFrame({"protein": df.loc[is_noninteraction, "bait"], "n_correct": n_correct_per_pair[is_noninteraction]})
            .groupby("protein")["n_correct"]
            .agg(n_correct_noninteraction="sum", n_pairs_noninteraction="count")
        )

        merged = inter_stats.join(non_stats, how="inner").reset_index()
        # denominator = n_pairs_for_protein * n_models_this_permutation; summed across permutations
        # downstream this becomes n_pairs * n_permutations * n_models_per_permutation = n_pairs * 60.
        merged["n_total_interaction"] = merged["n_pairs_interaction"] * n_models
        merged["n_total_noninteraction"] = merged["n_pairs_noninteraction"] * n_models
        merged["dataset"] = snakemake.wildcards.dataset
        merged["network_type"] = snakemake.wildcards.network_type
        merged["esm_model"] = snakemake.wildcards.esm_model
        merged["permutation"] = snakemake.wildcards.permutation
        merged["train_type"] = train_type
        frames.append(merged)

    result = pd.concat(frames, ignore_index=True)
    result = result[[
        "dataset", "network_type", "esm_model", "permutation", "train_type",
        "protein", "n_correct_interaction", "n_total_interaction",
        "n_correct_noninteraction", "n_total_noninteraction",
    ]]
    result.to_csv(snakemake.output.protein_accuracy, sep="\t", index=False)
