import numpy as np
import pandas as pd
from scipy.special import expit, logit


def main():
    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"

    pair_directional_df = pd.read_parquet(
        snakemake.input.pair_posterior,
        columns=[prot_a, prot_b, "pubmed_id", "log_odds", "llr_detected", "llr_not_detected", "n_observed", "n_tested"]
        )
    global_fit_df = pd.read_csv(snakemake.input.global_fit, sep="\t")

    pi = global_fit_df.loc[0, "pi"].astype(float)
    a, b = pair_directional_df[prot_a].values, pair_directional_df[prot_b].values
    pair_directional_df[prot_a] = np.where(a <= b, a, b)
    pair_directional_df[prot_b] = np.where(a <= b, b, a)

    pair_directional_df = pair_directional_df.groupby([prot_a, prot_b]).agg(
        log_odds=("log_odds", "sum"),
        llr_detected=("llr_detected", "sum"),
        llr_not_detected=("llr_not_detected", "sum"),
        n_observed=("n_observed", "sum"),
        n_tested=("n_tested", "sum"),
        pubmed_id=("pubmed_id", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),  # union of both directions' studies
        directions=("log_odds", "size"),
        ).reset_index()
    pair_directional_df.loc[pair_directional_df.directions > 1, "log_odds"] -= logit(pi) # remove the doubled prior
    pair_directional_df["posterior"] = expit(pair_directional_df.log_odds)

    pair_directional_df.to_parquet(snakemake.output.pair_posterior, index=False)


if __name__ == "__main__":
    main()
