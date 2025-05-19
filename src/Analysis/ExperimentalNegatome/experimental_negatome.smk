
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
        prior = 0.1
    input:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv"
    output:
        experimental_negatome = "work_folder/inferred_search_space/aggregated/methods/threshold_negatome.csv",
        hcl = "work_folder/inferred_search_space/aggregated/methods/high_confidence.csv"
    run:
        inferred_negative_df = pd.read(
            input.method_aggregate,
            sep="\t"
        )
        inferred_negative_df[
            (inferred_negative_df["n_observed"] == 0) &
            (inferred_negative_df["n_tested"] > params.min_observations)
        ].to_csv(output.experimental_negatome, sep="\t", index=False)

