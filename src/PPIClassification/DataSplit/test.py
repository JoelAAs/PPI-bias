import numpy as np
import pandas as pd

positive_data = "work_folder/per_gene/subsets/validation/ms_limit_0.15_maxpos_pos.csv"
negative_data = "work_folder/per_gene/subsets/validation/ms_limit_3_poslim_0.15_maxpos_neg.csv"

positive_df = pd.read_csv(positive_data, sep="\t")
negative_df = pd.read_csv(negative_data, sep="\t")

positive_df = positive_df.iloc[:, 0:2]
negative_df = negative_df.iloc[:, 0:2]

positive_df.columns = ["bait", "prey"]
negative_df.columns = ["bait", "prey"]

positive_df = positive_df[
    (positive_df["bait"].isin(negative_df["bait"])) |
    (positive_df["prey"].isin(negative_df["prey"]))
    ]  # remove negative edges with no negative df representation

negative_df = negative_df[
    (negative_df["bait"].isin(positive_df["bait"])) |
    (negative_df["prey"].isin(positive_df["prey"]))
    ]

all_targets = {f"{g}_prey" for g in positive_df["prey"].unique()} | {f"{g}_bait" for g in positive_df["bait"].unique()}
gene_idx = {gene: index + 1 for index, gene in enumerate(all_targets)}

# gene_idx[0] is not in list
target_list = np.zeros(len(gene_idx) + 1)
for gene, c in positive_df.groupby("prey", as_index=False).size().values:
    target_list[gene_idx[f"{gene}_prey"]] = c

for gene, c in positive_df.groupby("bait", as_index=False).size().values:
    target_list[gene_idx[f"{gene}_bait"]] = c

negative_df["prey_idx"] = negative_df["prey"].apply(lambda gene: gene_idx.get(f"{gene}_prey", 0))
negative_df["bait_idx"] = negative_df["bait"].apply(lambda gene: gene_idx.get(f"{gene}_bait", 0))

chosen_edges_mask = np.zeros(negative_df.shape[0])
bait_idx = negative_df['bait_idx'].values
prey_idx = negative_df['prey_idx'].values

for i in range(negative_df.shape[0]):
    print(f"{i}/{negative_df.shape[0]}")
    bait, prey = bait_idx[i], prey_idx[i]
    if target_list[bait] > 0 and target_list[prey] > 0:
        chosen_edges_mask[i] = 1
        target_list[bait] -= 1
        target_list[prey] -= 1

for i in range(negative_df.shape[0]):
    print(f"{i}/{negative_df.shape[0]}")
    if not chosen_edges_mask[i]:
        bait, prey = bait_idx[i], prey_idx[i]
        value_bait = target_list[bait]
        value_prey = target_list[prey]

        score_delta = (
                value_bait - (value_bait - 1 if bait != 0 else 0)
                + value_prey - (value_prey - 1 if prey != 0 else 0)
        )
        if

for i in remaining_edges:

while possible_decrease:
    for score_order in [1, 2]:
        remaining_edges_score = remaining_edges
        for i in remaining_edges:
            bait, prey = negative_df.iloc[i]
            val_bait = remaining_targets.get(f"{bait}_bait", 0)
            val_prey = remaining_targets.get(f"{prey}_prey", 0)
