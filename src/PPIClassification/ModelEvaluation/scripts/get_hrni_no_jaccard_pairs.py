import itertools
import re

import pandas as pd

MODEL_COLUMN_RE = re.compile(r"^train_(?P<train_type>hrni|no)_poslim(?P<pos_limit>all|[0-9.]+)_neglim(?P<neg_limit>[0-9.]+)$")


def jaccard(set_a, set_b):
    union = set_a | set_b
    if not union:
        return float("nan")
    return len(set_a & set_b) / len(union)


def predicted_positive_pairs(df, column):
    y_pred = df[column] >= 0.5
    return set(df.loc[y_pred, "pair_id"])


permutation_configs = {}
for pred_file in snakemake.input.predictions:
    permutation = pred_file.split("/permuted/")[1].split("/")[0]
    df = pd.read_csv(pred_file, sep="\t", dtype={"bait": "string", "prey": "string"})
    df = df[df["dataset"] == "Negative_HRNI"].copy()
    # network_type is undirectional: bait/prey is duplicated in both orientations, so identify
    # an edge by its sorted bait/prey pair rather than the directed (bait, prey) tuple.
    df["pair_id"] = df[["bait", "prey"]].min(axis=1) + ":" + df[["bait", "prey"]].max(axis=1)

    by_config = {}
    for column in df.columns:
        match = MODEL_COLUMN_RE.match(column)
        if match:
            key = (match["train_type"], match["pos_limit"], match["neg_limit"])
            by_config[key] = predicted_positive_pairs(df, column)
    permutation_configs[permutation] = by_config

configs = sorted({key[1:] for by_config in permutation_configs.values() for key in by_config})
permutations = sorted(permutation_configs)

rows = []
for pos_limit, neg_limit in configs:
    for permutation in permutations:
        hrni_set = permutation_configs[permutation][("hrni", pos_limit, neg_limit)]
        no_set = permutation_configs[permutation][("no", pos_limit, neg_limit)]
        rows.append({
            "dataset": snakemake.wildcards.dataset,
            "network_type": snakemake.wildcards.network_type,
            "pos_limit": pos_limit,
            "neg_limit": neg_limit,
            "comparison": "hrni_vs_no",
            "permutation_a": permutation,
            "permutation_b": permutation,
            "jaccard": jaccard(hrni_set, no_set),
        })

    for perm_a, perm_b in itertools.combinations(permutations, 2):
        set_a = permutation_configs[perm_a][("no", pos_limit, neg_limit)]
        set_b = permutation_configs[perm_b][("no", pos_limit, neg_limit)]
        rows.append({
            "dataset": snakemake.wildcards.dataset,
            "network_type": snakemake.wildcards.network_type,
            "pos_limit": pos_limit,
            "neg_limit": neg_limit,
            "comparison": "no_vs_no",
            "permutation_a": perm_a,
            "permutation_b": perm_b,
            "jaccard": jaccard(set_a, set_b),
        })

pd.DataFrame(rows).to_csv(snakemake.output.jaccard, sep="\t", index=False)