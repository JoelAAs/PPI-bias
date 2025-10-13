from scipy.stats import beta

checkpoint all_methods_filter_out:
    """
    Get upper and lower bound probability of detection given test/observations of each bait-prey combination
    """
    params:
        pseudo_n=config["pseudo_n"],
    input:
        method_aggregate=f"work_folder/{pn}/inferred_search_space/aggregated/methods/{{data}}_experimental_wise.csv"
    output:
        full_detection=f"work_folder/{pn}/analysis/POD/POD_{{data}}.csv"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df["gene_name_bait"] != inferred_negative_df["gene_name_prey"]
            ].copy()  # should be fixed later
        global_pod = inferred_negative_df["n_observed"].sum() / inferred_negative_df[
            "n_tested"].sum()  # if a test is made, probability of interaction
        prior_alpha = params.pseudo_n * global_pod
        prior_beta = params.pseudo_n - prior_alpha

        inferred_negative_df["alpha_post"] = prior_alpha + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"] = prior_beta + inferred_negative_df["n_tested"] - inferred_negative_df[
            "n_observed"]

        inferred_negative_df["p"] = inferred_negative_df["alpha_post"] / (
                inferred_negative_df["alpha_post"] + inferred_negative_df["beta_post"])

        inferred_negative_df["lower_bound_pod"] = beta.ppf(0.025,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["upper_bound_pod"] = beta.ppf(0.975,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])

        inferred_negative_df.to_csv(
            output.full_detection,
            sep="\t",
            index=False
        )


rule differential_detected_flat_negatome:
    #  TODO: Evaluate if this is used or useful
    """
    Flat interaction/non-interactions given cl specific prey-detection 
    """
    input:
        differential_interactions_filtered="work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_filtered.csv",
        experimental_negatome="work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hci="work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv"
    output:
        cl_negatome_cell="work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/threshold_negatome.csv",
        cl_hci="work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"
    run:
        df_diff = pd.read_csv(input.differential_interactions_filtered,sep="\t")
        df_nega = pd.read_csv(input.experimental_negatome,sep="\t")
        df_hci = pd.read_csv(input.hci,sep="\t")

        df_nega = df_nega.merge(
            df_diff,
            on="gene_name_prey"
        ).to_csv(output.cl_negatome_cell,sep="\t")
        df_hci = df_hci.merge(
            df_diff,
            on="gene_name_prey"
        ).to_csv(output.cl_hci,sep="\t")
