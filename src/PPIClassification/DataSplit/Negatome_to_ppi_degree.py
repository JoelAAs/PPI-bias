import datetime
import numpy as np
import pandas as pd
import argparse

def draw_and_update(bp_matrix, row_target, column_target, n_draws):
    all_row_index = np.zeros(n_draws)
    all_column_index = np.zeros(n_draws)
    all_row_errors = np.zeros(n_draws)
    all_column_errors = np.zeros(n_draws)

    for i in range(n_draws):
        row_sum = bp_matrix.sum(axis=1)
        column_sum = bp_matrix.sum(axis=0)
        total_sum = column_sum.sum()
        row_probability = row_sum - row_target * total_sum
        column_probability = column_sum - column_target * total_sum
        row_probability[row_probability < 0] = 0
        column_probability[column_probability < 0] = 0

        if sum(row_probability) == 0 or sum(column_probability) == 0:
            raise ValueError("No possible bait/prey choices")

        prob_matrix = row_probability[:, None] * column_probability[None, :]
        possible_choice_matrix = bp_matrix * prob_matrix
        flat_prob = possible_choice_matrix.ravel()
        flat_prob = flat_prob / flat_prob.sum()
        cum_probs = np.cumsum(flat_prob / flat_prob.sum())
        r = np.random.rand()
        flat_idx = np.searchsorted(cum_probs, r, side='right')

        # flat_idx = np.random.choice(range(len(flat_prob)), p=flat_prob / flat_prob.sum()) # too slow
        row_idx = int(flat_idx/possible_choice_matrix.shape[1])
        col_idx = flat_idx - row_idx*possible_choice_matrix.shape[1]
        row_sum[row_idx] -= 1
        column_sum[col_idx] -= 1

        delta_row = row_sum - row_target * total_sum
        delta_column = column_sum - column_target * total_sum
        row_error = np.abs(delta_row).sum()
        col_error = np.abs(delta_column).sum()
        bp_matrix[row_idx, col_idx] = 0

        all_row_index[i] = row_idx
        all_column_index[i] = col_idx
        all_row_errors[i] = row_error
        all_column_errors[i] = col_error

    return all_row_index, all_column_index, all_row_errors, all_column_errors, bp_matrix


def subset_negative_set(negative_bait_prey_df, positive_bait_prey_df):
    baits = set(negative_bait_prey_df["bait"]) | set(positive_bait_prey_df["bait"])
    bait_idx = {bait: i for i, bait in enumerate(baits)}
    idx_bait = {value: key for key, value in bait_idx.items()}

    all_prey = set(negative_bait_prey_df["prey"]) | set(positive_bait_prey_df["prey"])
    prey_idx = {prey: i for i, prey in enumerate(all_prey)}
    idx_prey = {value: key for key, value in prey_idx.items()}

    pos_bp_matrix = np.zeros((len(baits), len(all_prey)))
    neg_bp_matrix = np.zeros((len(baits), len(all_prey)))

    for bp_matrix, edge_df in zip(
            [neg_bp_matrix, pos_bp_matrix],
            [negative_bait_prey_df, positive_bait_prey_df]):
        for b, p in edge_df.values:
            bp_matrix[bait_idx[b], prey_idx[p]] = 1

    target_row_frequency = pos_bp_matrix.sum(axis=1) / pos_bp_matrix.sum()
    target_col_frequency = pos_bp_matrix.sum(axis=0) / pos_bp_matrix.sum()
    print("starting picking")
    with open("test.csv", "w") as w:
        w.write(f"bait\tprey\trow_error\tcol_error\n")
        n = 1000
        row_error = 10000
        col_error = 10000
        while (row_error + col_error) > 10000:
            s = datetime.datetime.now()
            batch_bait_idx_dropped, batch_prey_idx_dropped, batch_row_error, batch_col_error, neg_bp_matrix = draw_and_update(
                neg_bp_matrix, target_row_frequency, target_col_frequency, n)
            for bait_idx_dropped, prey_idx_dropped, row_error, col_error in zip(*
                    [batch_bait_idx_dropped, batch_prey_idx_dropped, batch_row_error, batch_col_error]):
                w.write(f"{idx_bait[bait_idx_dropped]}\t{idx_prey[prey_idx_dropped]}\t{row_error}\t{col_error}\n")
            e = datetime.datetime.now()
            print(f"{(e-s).seconds} seconds per {n} samples")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data
    model_name = args.model_name
    output_csv = args.embedding_csv



    df_ms = pd.read_csv("work_folder/per_gene/analysis/POD/POD_ms.csv", sep="\t")
    df_neg = df_ms[(df_ms["n_tested"] > 3) & (df_ms["n_observed"] == 0)]
    df_neg = df_neg[["gene_name_bait", "gene_name_prey"]].copy()
    df_neg.columns = ["bait", "prey"]

    df_pos = df_ms[df_ms["lower_bound_pod"] > 0.2]
    df_pos = df_pos[["gene_name_bait", "gene_name_prey"]].copy()
    df_pos.columns = ["bait", "prey"]

    subset_negative_set(df_neg, df_pos)