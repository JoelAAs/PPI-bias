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



rule all_methods_filter_out:
    params:
        min_observations = 4,
        pseudo_n = 1,
        HCL_fraction = 0.8
    input:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv"
    output:
        experimental_negatome = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hcl = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv"
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
        inferred_negative_df[inferred_negative_df["p"] > params.HCL_fraction].to_csv(
            output.hcl,
            sep="\t",
            index=False
        )
        inferred_negative_df[
            inferred_negative_df["p"] < prior_alpha/(prior_alpha + prior_beta + params.min_observations)
        ].to_csv(output.experimental_negatome, sep="\t", index=False)


rule create_cell_line_negatome_HCL:
    input:
        differential_interactions_filtered = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_filtered.csv",
        experimental_negatome= "work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hcl= "work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv"
    output:
        cl_negatome_cell = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/threshold_negatome.csv",
        cl_hcl = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"
    run:
        df_diff = pd.read_csv(input.differential_interactions_filtered, sep="\t")
        df_nega = pd.read_csv(input.experimental_negatome, sep="\t")
        df_hcl = pd.read_csv(input.hcl, sep="\t")


        df_nega = df_nega.merge(
            df_diff,
            by="gene_name_prey"
        ).to_csv(output.cl_negatome, sep="\t")
        df_hcl = df_hcl.merge(
            df_diff,
            by="gene_name_prey"
        ).to_csv(output.cl_hcl, sep="\t")