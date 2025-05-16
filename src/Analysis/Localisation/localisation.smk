import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


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
        localisation = "work_folder/inferred_search_space/analysis/localisation/logistic_method.csv"
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

        full_df["success_ratio"] = full_df["n_observed"]/full_df["n_tested"]
        formula = 'success_ratio ~ localisation_match + method'
        p = []

        localisations = full_df["localisation_bait"].unique()
        for localisation in localisations:
            loc_ss = full_df[full_df["localisation_bait"] == localisation]
            model = smf.glm(
                formula=formula, data=loc_ss,
                family=sm.families.Binomial(),
                freq_weights=loc_ss['n_tested']).fit()
            results = pd.concat([model.params, model.pvalues])
            results.index = [
                "intercept_coef", "match_coef", "method_coef",
                "intercept_pvalue", "match_pvalue", "method_pvalue"
            ]
            p.append(results)

        df_results = pd.DataFrame(p)
        df_results["localisation"] = localisations
        df_results.to_csv(output.localisation, sep="\t", index=False)