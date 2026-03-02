
def reformat_long_cl(df, pattern, col_name, id_column="gene_name_prey"):
    df_col = [c for c in df.columns if pattern in c]
    df_long = pd.melt(
        df[[id_column,] + df_col],
        id_vars=id_column,var_name="CVCL",value_name=col_name
    )
    df_long["CVCL"] = df_long["CVCL"].apply(lambda x: x.replace(pattern,""))

    return df_long

def nested_dict():
    return defaultdict(nested_dict)

def get_input_for_aggregation(wc, filename):
    CL_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_cell_line").output[0]

    cl_df = pd.read_csv(filename, sep="\t")
    cl_df = cl_df[
        cl_df[f"gene_name_bait"] != cl_df[f"gene_name_prey"]
        ] # remove bait-bait

    cl_df = cl_df[
        ~cl_df[[
            "gene_name_bait", "gene_name_prey",
            "pubmed_id", "detection_method", "CVCL"
        ]].duplicated(keep="first")] # remove isoforms

    cl_df = cl_df[["pubmed_id", "detection_method", "CVCL"]]
    cl_df = cl_df[cl_df.duplicated(keep=False)]

    expected_input = {
        f"{CL_FOLDER}/{pubmed_id}_{detection_method}_{cl_id}.csv"
        for _, (pubmed_id, detection_method, cl_id) in cl_df.iterrows()
    }
    return expected_input


rule aggregate_inferred_studies_cell_line:
    """
    Aggregate experiments assuming that any prey observed in studies is tested against all baits 
    """
    input:
        ppi_file = f"work_folder{pn}/formated/bait_prey_CVCL.csv",
        cl_pids = lambda wc: get_input_for_aggregation(wc, config["cell_line_ppis"])
    output:
        cell_line_counts = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    run:
        aggregate_inferred_experiments(input.cl_pids, output.cell_line_counts, cl=True)

rule infer_bait_wise_tests_cell_line:
    """
    Expand tests under the assumption that any prey that has been seen is tested in all 
    other experiments with the same bait.
    Params:
        remove_single_ppi_papers = True # Remove single ppi studies, assuming that they only focused on one interaction
    """
    params:
        remove_single_ppi_papers = config["remove_single_publications"]
    input:
        df = f"work_folder{pn}/formated/bait_prey_CVCL.csv",

    output:
        baitwise_infered = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_bait_wise.csv"
    run:
        ppi_df = pd.read_csv(input.df, sep="\t")

        
        # TODO: remove isoforms and bait-bait interactions
        ppi_df = ppi_df[
            ~ppi_df[["gene_name_bait", "gene_name_prey", "pubmed_id", "detection_method", "CVCL"]].duplicated()
        ].copy()
        ppi_df["study_id"] = ppi_df.apply(lambda row: str(row["pubmed_id"]) + row["detection_method"], axis=1)

        if params.remove_single_ppi_papers:
            single_experiment = ppi_df.groupby("study_id",as_index=False).size()
            single_experiment = single_experiment[single_experiment["size"] == 1]["study_id"]
            ppi_df = ppi_df[~ppi_df["study_id"].isin(single_experiment)]
        n_experiments = ppi_df.groupby(["gene_name_bait", "CVCL"], as_index=False)["study_id"].nunique()
        n_experiments = n_experiments.rename({
            "study_id": "n_tested"
        }, axis=1)

        combinations = ppi_df[~ppi_df[["gene_name_bait", "gene_name_prey"]].duplicated()][["gene_name_bait", "gene_name_prey"]]
        combinations = combinations.merge(n_experiments, on="gene_name_bait")
        n_observed = ppi_df.groupby(["gene_name_bait", "gene_name_prey", "CVCL"], as_index=False).size()
        n_observed = n_observed.rename({
            "size": "n_observed"
        }, axis=1)
        full_df = combinations.merge(n_observed, on=["gene_name_bait", "gene_name_prey", "CVCL"], how="left").fillna(0)

        tested = full_df[[
            'gene_name_bait', 'gene_name_prey', 'cl_id', 'n_tested'
        ]].pivot_table(
            index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',
            values='n_tested',
            fill_value=0)
        tested["total"] = tested.sum(axis=1)
        tested = tested.rename({
            column: f"{column}_tested" for column in tested.columns
        },axis=1).reset_index()

        observed = full_df[[
            'gene_name_bait', 'gene_name_prey', 'cl_id', 'n_observed'
        ]].pivot_table(
            index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',
            values='n_observed',
            fill_value=0)
        observed["total"] = observed.sum(axis=1)
        observed = observed.rename({
            column: f"{column}_observed" for column in observed.columns
        },axis=1).reset_index()

        tested_observed_wide = observed.merge(tested, on = ['gene_name_bait', 'gene_name_prey'])
        tested_observed_wide.to_csv(output.baitwise_infered, sep="\t", index=False)


rule test_prey_probability:
    """
    Assuming the probability of observing an interaction is the same regardless of experiment.
    As P(p|b1) = P(p|b2) if P(p|b1), P(p|b2) != 0 
    We test selected cell lines on a prey basis against all other cell lines (not only selected)  
    """
    params:
        selected_celllines = config["selected_cell_lines"]
    input:
        bait_wise_inferred = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_bait_wise.csv"
    output:
        prey_bait_wise_tested = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_tested.csv"
    run:
        bait_wise_inferred      = pd.read_csv(input.bait_wise_inferred, sep="\t")
        bait_wise_inferred_prey = bait_wise_inferred.groupby("gene_name_prey", as_index=False).sum()
        columns_to_keep = ["gene_name_prey",] + [
            f"{selected_cl}_observed" for selected_cl in params.selected_celllines
        ] + [
            f"{selected_cl}_tested" for selected_cl in params.selected_celllines
        ] + ["total_tested", "total_observed"]
        bait_wise_inferred_prey = bait_wise_inferred_prey[columns_to_keep]

        def _fisher_exact(row, c_cl):
            c_observed = row[f"{c_cl}_observed"]
            c_tested = row[f"{c_cl}_tested"]
            c_non_observed = c_tested - c_observed
            other_observed = row["total_observed"] - c_observed
            other_tested = row["total_tested"] - c_tested
            other_non_observed = other_tested - other_observed

            c_table = [
                [c_observed, c_non_observed],
                [other_observed, other_non_observed]
            ]

            # if c_tested == 0 | other_observed == 0:
            #     OR = np.nan
            OR, p_value  = fisher_exact(c_table)
            return p_value, OR

        for cl in params.selected_celllines:
            bait_wise_inferred_prey[f"{cl}_p_value"], bait_wise_inferred_prey[f"{cl}_odds_ratio"] = zip(*bait_wise_inferred_prey.apply(
        lambda row: _fisher_exact(row, cl), axis=1))


        bait_wise_inferred_prey.to_csv(output.prey_bait_wise_tested , sep="\t", index=False)

rule filter_tests:
    """
    Don't got changing these parameters after specific results are found.
    That's cheating
    """
    params:
        min_total_tests = config["min_total_tests"],
        min_total_observed = config["min_total_observed"],
        part_or_difference_cutoff = config["part_or_difference_cutoff"],
        min_fdr_pval = 0.05
    input:
        prey_bait_wise_tested = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_tested.csv",
    output:
        differential_interactions_filtered = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_filtered.csv",
        differential_interactions_ploting  = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv"
    run:
        results_df = pd.read_csv(input.prey_bait_wise_tested, sep="\t")
        results_df = results_df[
            (results_df["total_observed"] > params.min_total_observed) &
            (results_df["total_tested"] > params.min_total_tests)
        ]
        p_columns = [c for c in results_df.columns if "_p_value" in c]
        p_value_df_long = pd.melt(
            results_df[["gene_name_prey"] + p_columns],
            id_vars='gene_name_prey',
            var_name="CVCL",
            value_name="p_value"
        )
        p_value_df_long["CVCL"] = p_value_df_long["CVCL"].apply(lambda x: x.replace("_p_value", ""))
        p_value_df_long["p_value_adjusted"] = false_discovery_control(p_value_df_long["p_value"], method="by")

        or_columns = [c for c in results_df.columns if "_odds_ratio" in c]
        or_df_long = pd.melt(
            results_df[["gene_name_prey"] + or_columns],
            id_vars='gene_name_prey',
            var_name="CVCL",
            value_name="odds_ratio"
        )
        or_df_long["CVCL"] = or_df_long["CVCL"].apply(lambda x: x.replace("_odds_ratio", ""))

        results_filtered_long = or_df_long.merge(p_value_df_long, on=["gene_name_prey", "CVCL"])
        results_filtered_long = results_filtered_long[
            (results_filtered_long["p_value_adjusted"] < params.min_fdr_pval)
        ]
        results_filtered_long.to_csv(output.differential_interactions_ploting, sep="\t", index=False)
        results_filtered_long = results_filtered_long[
                    (results_filtered_long["odds_ratio"] < 1/params.part_or_difference_cutoff) |
                    (results_filtered_long["odds_ratio"] > params.part_or_difference_cutoff)
        ]
        results_filtered_long.to_csv(output.differential_interactions_filtered, sep="\t", index=False)


rule create_cell_line_negatome_HCL:
    params:
        min_observations = 1,
        pseudo_n = 1,
        top_fraction = 0.6 # There isn't that many experiments when dividing on N_cl
    input:
        differential_interactions_filtered = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_filtered.csv",
        experiment_wise = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    output:
        cell_line_negatome = "work_folder/inferred_search_space/analysis/cell_line/negatome_cl.csv",
        hcl = "work_folder/inferred_search_space/analysis/cell_line/HCL_cl.csv"
    run:
        cl_diff = pd.read_csv(input.differential_interactions_filtered, sep = "\t")
        experimental_df = pd.read_csv(input.experiment_wise, sep = "\t")


        min_test_df = experimental_df[
            experimental_df["n_tested"] > params.min_observations
            ].copy()
        min_test_df["ratio"] = min_test_df["n_observed"] / min_test_df["n_tested"]
        mean_p = min_test_df["ratio"].mean()
        prior_alpha = params.pseudo_n * mean_p
        prior_beta = params.pseudo_n - prior_alpha

        experimental_df["alpha_post"] = prior_alpha + experimental_df["n_observed"]
        experimental_df["beta_post"] = prior_beta + experimental_df["n_tested"] - experimental_df["n_observed"]

        experimental_df["p"] = experimental_df["alpha_post"] / (
                experimental_df["alpha_post"] + experimental_df["beta_post"])

        hcl_cell_line = experimental_df[experimental_df["p"] > params.top_fraction]
        hcl_cell_line = hcl_cell_line.merge(
            cl_diff,
            on=["gene_name_prey", "CVCL"]
        )
        hcl_cell_line.to_csv(
            output.hcl,
            sep="\t",
            index=False
        )
        negatome = experimental_df[
            experimental_df["p"] < prior_alpha / (prior_alpha + prior_beta + params.min_observations)
            ]
        negatome = negatome.merge(
            cl_diff,
            on=["gene_name_prey", "CVCL"]
        )
        negatome.to_csv(output.cell_line_negatome,sep="\t",index=False)


rule marginalised_prey_probability:
    params:
        prior_strength = 1/3,
        selected_celllines = config["selected_cell_lines"],
    input:
        bait_wise_inferred = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_bait_wise.csv"
    output:
        bait_based_prior = "work_folder/inferred_search_space/analysis/cell_line/bait_prior.csv",
        bait_based_prior_long= "work_folder/inferred_search_space/analysis/cell_line/bait_prior_long.csv",
    run:
        df_tests = pd.read_csv(input.bait_wise_inferred, sep="\t")
        df_tests["total_not_observed"] = df_tests[f"total_tested"] - df_tests[f"total_observed"]
        df_tests["prior_alpha"] = params.prior_strength * df_tests[f"total_observed"]
        df_tests["prior_beta"] = params.prior_strength * df_tests["total_not_observed"]

        for cell_line in params.selected_celllines:
            df_tests[f"{cell_line}_post_alpha"] = df_tests["prior_alpha"] + df_tests[f"{cell_line}_observed"]
            df_tests[f"{cell_line}_post_beta"]  = df_tests["prior_beta"]  + df_tests[f"{cell_line}_tested"] - df_tests[f"{cell_line}_observed"]
            df_tests[f"p_bait_{cell_line}"] = df_tests[f"{cell_line}_post_alpha"]/(df_tests[f"{cell_line}_post_beta"] + df_tests[f"{cell_line}_post_alpha"])
            df_tests[f"{cell_line}_bait_weight"] = df_tests[f"{cell_line}_tested"]/df_tests.groupby("gene_name_prey")[f"{cell_line}_tested"].transform("sum")
            df_tests[f"{cell_line}_p_bait_weighted"]  = df_tests[f"p_bait_{cell_line}"] * df_tests[f"{cell_line}_bait_weight"] # Weighted mean of P

        cl_alpha_prior_cols = [
            f"{cell_line}_observed" for cell_line in params.selected_celllines
        ] + [
            f"{cell_line}_tested" for cell_line in params.selected_celllines
        ] + [
            f"{cell_line}_p_bait_weighted" for cell_line in params.selected_celllines
        ] + [
            "total_observed",
            "total_tested"
        ]
        df_prey_probability = df_tests.groupby("gene_name_prey", as_index=False)[cl_alpha_prior_cols].sum()
        df_prey_probability.to_csv(
            output.bait_based_prior,
            sep="\t",
            index=False
        )

        n_tests = df_prey_probability[[
            "gene_name_prey",
            "total_tested",
            "total_observed"
        ]].copy()
        del df_prey_probability["total_tested"]
        del df_prey_probability["total_observed"]

        p_long = reformat_long_cl(df_prey_probability, "_p_bait_weighted", "prey_p")
        tested_long = reformat_long_cl(df_prey_probability,"_tested","n_tested")
        observed_long = reformat_long_cl(df_prey_probability,"_observed","n_observed")

        df_long = p_long.merge(
            tested_long.merge(
                observed_long.merge(
                    n_tests, on="gene_name_prey"
                ), on = ["gene_name_prey", "CVCL"]
            ),on=["gene_name_prey", "CVCL"]
        )

        df_long.to_csv(
            output.bait_based_prior_long,
            sep="\t",
            index=False
        )



