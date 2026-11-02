import pandas as pd
from scipy.special import logit


def manin():
    pair_directional_df = pd.read_parquet(
        snakemake.input.pair_posterior,
        columns=[prot_a, prot_b, "log_odds", "n_observed", "n_tested"]
        )
    global_fit_df = pd.read_csv(snakemake.input.global_fit, sep="\t")
    
    pi = global_fit_df.loc[0, "pi"].astype(float)
    pair_directional_df["prot_a"] = np.minimum(pair_directional_df[prot_a], pair_directional_df[prot_b])
    pair_directional_df["prot_b"] = np.maximum(pair_directional_df[prot_a], pair_directional_df[prot_b])
    pair_directional_df[prot_a] = pair_directional_df["prot_a"]
    pair_directional_df[prot_b] = pair_directional_df["prot_b"]
    del pair_directional_df["prot_a"]; del pair_directional_df["prot_b"]
    
    pair_directional_df = pair_directional_df.groupby([prot_a, prot_b]).agg(
        {"log_odds": "sum",
         "n_observed": "sum",
         "n_tested": "sum",
         "directions": "count"}).reset_index()
    pair_directional_df[pair_directional_df.directions > 1]["log_odds"] -= logit(pi) # remove the doubled prior
    
    pair_directional_df.to_parquet(snakemake.output.pair_posterior, index=False)
    
    