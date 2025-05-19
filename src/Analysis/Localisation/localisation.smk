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
        left_on="gene_name_bait",
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
        direction = "->".join(loc_ss[test_column])
        p.append({
            "OR": or_ln,
            "p_value": p_value,
            "testdirection": direction
        })

    df_results = pd.DataFrame(p)
    df_results["localisation"] = localisations
    return  df_results



rule method_comparison:
    """
    Compare localisation
    """
    params:
        localisation_csv = config["localisation_file"],
        min_localisation_genes = 100
    input:
        multi_method_ms  = "work_folder/inferred_search_space/aggregated/multi_methods/ms_experimental_wise.csv",
        multi_method_y2h = "work_folder/inferred_search_space/aggregated/multi_methods/y2h_experimental_wise.csv"
    output:
        method_diff_localisation = "work_folder/inferred_search_space/analysis/localisation/same_localisation_method_diff.csv",
        ms_diff_localisation = "work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
        y2h_diff_localisation= "work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv"
    run:
        df_ms  = pd.read_csv(input.multi_method_ms,  sep="\t")
        df_y2h = pd.read_csv(input.multi_method_y2h, sep="\t")

        df_localisation = pd.read_csv(params.localisation_csv, sep="\t")
        localisation_count = df_localisation.groupby("localisation", as_index=False).size()
        keep_localisations = localisation_count[
            localisation_count["size"] > params.min_localisation_genes
        ]["localisation"]
        df_localisation = df_localisation[
            df_localisation["localisation"].isin(keep_localisations)
        ]

        df_ms = add_localisation(df_ms, df_localisation)
        df_ms = df_ms.groupby(['localisation_bait', 'localisation_match'], as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_ms["method"] = "ms"

        df_y2h = add_localisation(df_y2h, df_localisation)
        df_y2h = df_y2h.groupby(['localisation_bait', 'localisation_match'], as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_y2h["method"] = "y2h"

        full_df = pd.concat([df_y2h, df_ms])
        unseen  = full_df.groupby("localisation_bait", as_index=False).size()
        localisation_unseen = unseen[unseen["size"] != 4]["localisation_bait"]
        full_df = full_df[~full_df["localisation_bait"].isin(localisation_unseen)]
        full_df["not_observed"] = full_df["n_tested"] - full_df["n_observed"]
        localisations = full_df["localisation_bait"].unique()

        diff_true = fisher_exact_to_df(
            full_df[full_df["localisation_matched"] == True],
            localisations,
            "method"
        )
        diff_true.to_csv(output.method_diff_localisation, sep="\t", index=False)

        diff_ms = fisher_exact_to_df(
            full_df[full_df["method"] == "ms"],
            localisations,
            "localisation_matched"
        )
        diff_ms.to_csv(output.ms_diff_localisation, sep="\t", index=False)

        diff_y2h = fisher_exact_to_df(
            full_df[full_df["method"] == "y2h"],
            localisations,
            "localisation_matched"
        )
        diff_y2h.to_csv(output.y2h_diff_localisation, sep="\t",index=False)
