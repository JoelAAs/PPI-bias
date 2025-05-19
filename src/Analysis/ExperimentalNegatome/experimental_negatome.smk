from openpyxl.styles.builtins import output

rule aggregate_pids:
    """
    Aggregate data form studies of the same method
    """
    input:
        method = expand(
            "work_folder/inferred_search_space/aggregated/multi_methods/{multi_method}_experimental_wise.csv",
            multi_method=["ms", "y2h"]
        )
    output:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv"
    run:
        aggregate_inferred_experiments(input.method, output.method_aggregate)



rule filter_out:
    params:
        min_observations = 4,
        pseudo_n = 1
    input:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv"
    output:
        experimental_negatome = "work_folder/inferred_search_space/aggregated/methods/threshold_negatome.csv",
        hcl = "work_folder/inferred_search_space/aggregated/methods/high_confidence.csv"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        min_test_df = inferred_negative_df[
            inferred_negative_df["n_tested"] > params.min_observations
        ]
        min_test_df["ratio"] = min_test_df["n_observed"]/min_test_df["n_tested"]
        mean_p = min_test_df["ratio"].mean()
        prior_alpha = params.pseudo_n*mean_p
        prior_beta  = params.pseudo_n - prior_alpha

        inferred_negative_df["alpha_post"] = prior_alpha + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"]  = prior_beta + min_test_df["n_tested"] - inferred_negative_df["n_observed"]

        inferred_negative_df["p"] = inferred_negative_df["alpha_post"]/(
                inferred_negative_df["alpha_post"] + inferred_negative_df["beta_post"])
        inferred_negative_df[inferred_negative_df["p"] > 0.90].to_csv(
            output.hcl,
            sep="\t",
            index=False
        )
        inferred_negative_df[
            inferred_negative_df["p"] < prior_alpha/(prior_alpha + prior_beta + params.min_observations)
        ].to_csv(output.experimental_negatome, sep="\t", index=False)
