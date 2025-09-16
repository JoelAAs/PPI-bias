import pandas as pd
from scipy.stats import fisher_exact


def add_localisation(ppi_df, localisation_df):
    ppi_df = ppi_df.merge(
        localisation_df,
        left_on="gene_name_bait",
        right_on="gene_name"
    )
    del ppi_df["gene_name"]
    ppi_df = ppi_df.merge(
        localisation_df,
        left_on="gene_name_prey",
        right_on="gene_name",
        suffixes=("_bait", "_prey")
    )
    del ppi_df["gene_name"]
    ppi_df["localisation_match"] = ppi_df["localisation_bait"] == ppi_df["localisation_prey"]

    return ppi_df

def fisher_exact_to_df(df_in, localisations, test_column):
    p = []
    for localisation in localisations:
        loc_ss         = df_in[df_in["localisation_bait"] == localisation]
        or_ln, p_value = fisher_exact(loc_ss[["n_observed", "not_observed"]])
        direction = "->".join(map(str,loc_ss[test_column]))
        p.append({
            "OR": or_ln,
            "p_value": p_value,
            "test_direction": direction
        })

    df_results = pd.DataFrame(p)
    df_results["localisation"] = localisations
    return  df_results

