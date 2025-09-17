from localisation_support import *
from mean_distance_support import get_cumulative_sum


def get_probabilities(df, value_column, greater=True, min_samples=50):
    bins = df[value_column].unique()
    bins.sort()

    values = df[value_column].values
    idx_val = values.argsort()
    if greater:
        idx_val = idx_val[::-1]
        bins = bins[::-1]
    values = values[idx_val]
    probabilities = df["match_probability"].to_numpy()[idx_val]
    matches = df["localisation_match"].to_numpy()[idx_val]

    if not greater:
        bins = -bins
        values = -values

    expected = 0
    observed = 0
    previous = 0
    i = 0
    j = 0
    rows = [{}] * len(bins)
    for threshold in bins:
        while (i < len(values) and threshold <= values[i]) or i < min_samples:
            i += 1
        if previous != i:
            expected += probabilities[previous:i].sum()
            observed += matches[previous:i].sum()
            previous = i
            rows[j] = {
                "limit": value_column,
                "value": (threshold if greater else -threshold),
                "expected": expected,
                "observed": observed,
                "number_of_pairs": i
            }
            j += 1
    return pd.DataFrame(rows).dropna()


rule method_comparison:
    """
    Compare localisation
    """
    params:
        localisation_csv=config["localisation_file"],
        min_localisation_genes=100
    input:
        multi_method_ms="work_folder/inferred_search_space/aggregated/multi_methods/ms_experimental_wise.csv",
        multi_method_y2h="work_folder/inferred_search_space/aggregated/multi_methods/y2h_experimental_wise.csv"
    output:
        method_diff_localisation="work_folder/inferred_search_space/analysis/localisation/same_localisation_method_diff.csv",
        ms_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
        y2h_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_y2h.csv"
    run:
        df_ms = pd.read_csv(input.multi_method_ms,sep="\t")
        df_y2h = pd.read_csv(input.multi_method_y2h,sep="\t")

        df_localisation = pd.read_csv(params.localisation_csv,sep="\t")
        localisation_count = df_localisation.groupby("localisation",as_index=False).size()
        keep_localisations = localisation_count[
            localisation_count["size"] > params.min_localisation_genes
            ]["localisation"]
        df_localisation = df_localisation[
            df_localisation["localisation"].isin(keep_localisations)
        ]

        df_ms = add_localisation(df_ms,df_localisation)
        df_ms = df_ms.groupby(['localisation_bait', 'localisation_match'],as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_ms["method"] = "ms"

        df_y2h = add_localisation(df_y2h,df_localisation)
        df_y2h = df_y2h.groupby(['localisation_bait', 'localisation_match'],as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_y2h["method"] = "y2h"

        full_df = pd.concat([df_y2h, df_ms])
        unseen = full_df.groupby("localisation_bait",as_index=False).size()
        localisation_unseen = unseen[unseen["size"] != 4]["localisation_bait"]
        full_df = full_df[~full_df["localisation_bait"].isin(localisation_unseen)]
        full_df["not_observed"] = full_df["n_tested"] - full_df["n_observed"]
        localisations = full_df["localisation_bait"].unique()

        diff_true = fisher_exact_to_df(
            full_df[full_df["localisation_match"] == True],
            localisations,
            "method"
        )
        diff_true.to_csv(output.method_diff_localisation,sep="\t",index=False)

        diff_ms = fisher_exact_to_df(
            full_df[
                full_df["method"] == "ms"
                ].sort_values('localisation_match',ascending=False),
            localisations,
            "localisation_match"
        )
        diff_ms.to_csv(output.ms_diff_localisation,sep="\t",index=False)

        diff_y2h = fisher_exact_to_df(
            full_df[
                full_df["method"] == "y2h"
                ].sort_values('localisation_match',ascending=False),
            localisations,
            "localisation_match"
        )
        diff_y2h.to_csv(output.y2h_diff_localisation,sep="\t",index=False)


rule accumulation_colocalisation:
    params:
        localisation_csv=config["localisation_file"]
    input:
        pod_data="work_folder/analysis/POD/POD_{data}.csv"
    output:
        localisation_annotated="work_folder/analysis/localisation/POD_{data}_localisation.csv",
        localisation_lesser="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_lesser.csv",
        localisation_greater="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_greater.csv"
    run:
        df_localisation = pd.read_csv(params.localisation_csv,sep="\t")
        bait_model = pd.read_csv(input.pod_data,sep="\t")
        proteins_tested = set(bait_model["gene_name_bait"].tolist() + bait_model["gene_name_prey"].tolist())
        df_localisation = df_localisation[df_localisation["gene_name"].isin(proteins_tested)]  # detectable
        n_possible = df_localisation.shape[0]
        random_match = pd.DataFrame([
            [loc, ((df_localisation["localisation"] == loc).sum() - 1) / n_possible] for loc in
            df_localisation["localisation"].unique()
        ],columns=("localisation_prey", "match_probability"))

        bait_model = add_localisation(bait_model,df_localisation)
        bait_model = bait_model.merge(random_match,on="localisation_prey")
        bait_model.to_csv(
            output.localisation_annotated,sep="\t",index=False)

        measurement_columns = ["match_probability", "localisation_match"]
        get_cumulative_sum(
            bait_model,
            value_column="lower_bound_pod",
            cumulative_columns=measurement_columns).to_csv(
            output.localisation_greater,sep="\t",index=False)
        get_cumulative_sum(
            bait_model,
            value_column="lower_bound_pod",
            cumulative_columns=measurement_columns,
            greater=False).to_csv(
            output.localisation_lesser,sep="\t",index=False)
