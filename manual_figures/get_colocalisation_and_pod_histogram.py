import math

import pandas as pd

data="y2h"

df_data = pd.read_csv(f"work_folder/per_gene/analysis/POD/POD_{data}.csv", sep="\t")
localisation_data = pd.read_csv(f"work_folder/per_gene/analysis/localisation/study_match_probability/expected/pairs_{data}_expected.csv", sep="\t")


df_data = df_data.merge(localisation_data, on="pair_id")
df_data_ss = df_data[df_data["n_observed"] != 0]

min_pod = df_data_ss["lower_bound_pod"].min()
expected = df_data_ss["match_probability"].sum()
observed = df_data_ss["localisation_match"].sum()

with open(f"exp_obs_{data}.txt") as w:
    w.write(f"{data}\n")
    w.write(f"Min pod: {min_pod}\n")
    w.write(f"expected: {expected}\n")
    w.write(f"Observed: {observed}\n")

logit = lambda x: math.log(x/(1-x))

df_data["upper_bound_pod_log"] = df_data["upper_bound_pod"].apply(logit)

pd.cut(df_data["upper_bound_pod_log"], 40).sort_values().to_csv(f"{data}_upper_log.csv", sep="\t")
pd.cut(df_data["lower_bound_pod"], 40).sort_values().to_csv(f"{data}_lower.csv", sep="\t")