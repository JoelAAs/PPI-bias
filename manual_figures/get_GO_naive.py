import math

import pandas as pd

data="y2h"

df_data = pd.read_csv(f"work_folder/per_gene/analysis/POD/POD_{data}.csv", sep="\t")
localisation_data = pd.read_csv(f"work_folder/per_gene/analysis/GO/pairs_{data}_jaccard.csv", sep="\t")


df_data = df_data.merge(localisation_data, on="pair_id")
df_data_ss = df_data[df_data["n_observed"] != 0]

min_pod = df_data_ss["lower_bound_pod"].min()

with open(f"go_{data}.txt") as w:
    w.write(f"{data}\n")
    w.write(f"Min pod: {min_pod}\n")
    w.write(f"Mean mf: {df_data_ss['sum_ji_mf'].sum() / df_data_ss['non_na_pairs_ji_mf'].sum()}\n")
    w.write(f"Mean bp: {df_data_ss['sum_ji_bp'].sum() / df_data_ss['non_na_pairs_ji_bp'].sum()}\n")
