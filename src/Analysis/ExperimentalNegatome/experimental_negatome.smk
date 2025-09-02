from scipy.stats import beta

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
    """
    TODO: OBS bait-bait still here
    """
    params:
        negatome_tested_threshold = config["negatome_tested_threshold"],
        pseudo_n = config["pseudo_n"],
        hci_fraction = config["hci_frac"]
    input:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv"
    output:
        experimental_negatome = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hci = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv",
        full_detection = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/p_estimated_protein_pairs.csv"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df["gene_name_bait"] != inferred_negative_df["gene_name_prey"]
        ].copy()  # should be fixed later
        global_pod = inferred_negative_df["n_observed"].sum()/inferred_negative_df["n_tested"].sum() # if a test is made, probability of interaction
        prior_alpha = params.pseudo_n*global_pod
        prior_beta  = params.pseudo_n - prior_alpha

        inferred_negative_df["alpha_post"] = prior_alpha + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"]  = prior_beta + inferred_negative_df["n_tested"] - inferred_negative_df["n_observed"]

        inferred_negative_df["p"] = inferred_negative_df["alpha_post"]/(
                inferred_negative_df["alpha_post"] + inferred_negative_df["beta_post"])

        inferred_negative_df["p_lower_ci"] = beta.ppf(0.025,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["p_upper_ci"] = beta.ppf(0.975,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])

        inferred_negative_df.to_csv(
            output.full_detection,
            sep="\t",
            index=False
        )

        inferred_negative_df[inferred_negative_df["p"] > params.hci_fraction].to_csv(
            output.hci,
            sep="\t",
            index=False
        )
        inferred_negative_df[
            inferred_negative_df["p"] < prior_alpha/(prior_alpha + prior_beta + params.negatome_tested_threshold)
        ].to_csv(output.experimental_negatome, sep="\t", index=False)


rule differential_detected_flat_negatome:
    """
    Flat interaction/non-interactions given cl specific prey-detection 
    """
    input:
        differential_interactions_filtered = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_filtered.csv",
        experimental_negatome= "work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hci= "work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv"
    output:
        cl_negatome_cell = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/threshold_negatome.csv",
        cl_hci = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"
    run:
        df_diff = pd.read_csv(input.differential_interactions_filtered, sep="\t")
        df_nega = pd.read_csv(input.experimental_negatome, sep="\t")
        df_hci  = pd.read_csv(input.hci, sep="\t")


        df_nega = df_nega.merge(
            df_diff,
            on="gene_name_prey"
        ).to_csv(output.cl_negatome_cell, sep="\t")
        df_hci = df_hci.merge(
            df_diff,
            on="gene_name_prey"
        ).to_csv(output.cl_hci, sep="\t")