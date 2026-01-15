import datetime
import numpy as np
import pandas as pd
import argparse


def draw_and_update(bp_matrix, row_target, column_target, subtractive, n_draws):
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
        if not subtractive:
            row_probability = -row_probability
            column_probability = -column_probability
        row_probability[row_probability < 0] = 0
        column_probability[column_probability < 0] = 0

        if sum(row_probability) == 0 or sum(column_probability) == 0:
            raise ValueError("No possible bait/prey choices")

        if bp_matrix.sum() == 0:
            raise ValueError("No PPIs left")

        prob_matrix = row_probability[:, None] * column_probability[None, :]
        possible_choice_matrix = bp_matrix * prob_matrix
        flat_prob = possible_choice_matrix.ravel()
        flat_prob = flat_prob / flat_prob.sum()
        cum_probs = np.cumsum(flat_prob / flat_prob.sum())
        r = np.random.rand()
        flat_idx = np.searchsorted(cum_probs, r, side='right')

        # flat_idx = np.random.choice(range(len(flat_prob)), p=flat_prob / flat_prob.sum()) # too slow
        row_idx = int(flat_idx / possible_choice_matrix.shape[1])
        col_idx = flat_idx - row_idx * possible_choice_matrix.shape[1]
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


def subset_negative_set(negative_bp_df, positive_bp_df, select_ppi_file, subtractive_bool, size_setting,
                        acceptable_error):
    baits = set(negative_bp_df["bait"]) & set(positive_bp_df["bait"])
    bait_idx = {bait: i for i, bait in enumerate(baits)}
    idx_bait = {value: key for key, value in bait_idx.items()}

    all_prey = set(negative_bp_df["prey"]) & set(positive_bp_df["prey"])
    prey_idx = {prey: i for i, prey in enumerate(all_prey)}
    idx_prey = {value: key for key, value in prey_idx.items()}

    pos_bp_matrix = np.zeros((len(baits), len(all_prey)), dtype=int)
    neg_bp_matrix = np.zeros((len(baits), len(all_prey)), dtype=int)

    # Remove all PPIs where bait and prey is not present in both neg/pos
    negative_bp_df = negative_bp_df[
        negative_bp_df["bait"].isin(baits) & negative_bp_df["prey"].isin(all_prey)]
    positive_bp_df = positive_bp_df[
        positive_bp_df["bait"].isin(baits) & positive_bp_df["prey"].isin(all_prey)]

    for bp_matrix, edge_df in zip(
            [neg_bp_matrix, pos_bp_matrix],
            [negative_bp_df, positive_bp_df]):
        for b, p in edge_df.values:
            bp_matrix[bait_idx[b], prey_idx[p]] = 1

    target_row_frequency = pos_bp_matrix.sum(axis=1) / pos_bp_matrix.sum()
    target_col_frequency = pos_bp_matrix.sum(axis=0) / pos_bp_matrix.sum()

    if size_setting == "equal":
        picked_limit = pos_bp_matrix.sum()
    else:
        picked_limit = neg_bp_matrix.sum()

    percent_degree_mean_error = round(pos_bp_matrix.sum() * acceptable_error)
    row_error = percent_degree_mean_error; col_error = percent_degree_mean_error
    print(f"Starting picking ppis with a aimed difference of at most: {percent_degree_mean_error}")
    n = 1000
    picked = 0
    with open(select_ppi_file, "w") as w:
        w.write(f"bait\tprey\trow_error\tcol_error\n")
        while (row_error + col_error) > percent_degree_mean_error and picked < picked_limit:
            s = datetime.datetime.now()
            batch_bait_idx_dropped, batch_prey_idx_dropped, batch_row_error, batch_col_error, neg_bp_matrix = draw_and_update(
                neg_bp_matrix, target_row_frequency, target_col_frequency, subtractive_bool, n)
            for bait_idx_dropped, prey_idx_dropped, row_error, col_error in zip(
                    *[batch_bait_idx_dropped,
                      batch_prey_idx_dropped,
                      batch_row_error,
                      batch_col_error]):
                w.write(f"{idx_bait[bait_idx_dropped]}\t{idx_prey[prey_idx_dropped]}\t{row_error}\t{col_error}\n")
            e = datetime.datetime.now()
            picked += n
            print(f"{(e - s).seconds} seconds per {n} samples")

    selected_ppi_df = pd.read_csv(select_ppi_file, sep="\t")
    selected_ppi_df["mean_error"] = selected_ppi_df["row_error"] + selected_ppi_df["col_error"]

    if size_setting == "equal" and not subtractive_bool:
        return selected_ppi_df.iloc[:picked_limit][["bait", "prey"]], positive_bp_df
    else:
        last_row = selected_ppi_df[selected_ppi_df["mean_error"] < percent_degree_mean_error]
        last_row_idx = last_row.index.tolist()[0]
        selected_ppi_df = selected_ppi_df.iloc[:last_row_idx]
        remove_ids = selected_ppi_df[["bait", "prey"]].apply(lambda x: ":".join(x)).tolist()
        negative_bp_df["id"] = negative_bp_df[["bait", "prey"]].apply(lambda x: ":".join(x))
        negative_bp_df = negative_bp_df[negative_bp_df["id"].isin(remove_ids)]
        return negative_bp_df[["bait", "prey"]], positive_bp_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--selected_ppis", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    parser.add_argument("--subtractive", default=False, help="Path to output csv file")
    parser.add_argument("--size", default="max", help="Path to output csv file")
    parser.add_argument("--accepted_error", default=0.1, required=True, help="Path to output csv file")

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data

    selected_ppi_file = args.selected_ppis

    balanced_negative = args.balanced_negative
    balanced_positive = args.balanced_positive
    subtractive = args.subtractive
    size = args.size
    accepted_error = args.accepted_error

    negative_bait_prey_df = pd.read_csv(negative_data, sep="\t")
    positive_bait_prey_df = pd.read_csv(positive_data, sep="\t")

    negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
    positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

    negative_bait_prey_df.columns = ["bait", "prey"]
    positive_bait_prey_df.columns = ["bait", "prey"]

    balanced_negative_df, balanced_positive_df =subset_negative_set(
        negative_bait_prey_df,
        positive_bait_prey_df,
        selected_ppi_file,
        subtractive,
        size,
        accepted_error
    )
    balanced_negative_df.to_csv(balanced_negative, sep="\t", index=False)
    balanced_positive_df.to_csv(balanced_positive, sep="\t", index=False)
