from scipy.stats import beta
import datetime

rule all_methods_filter_out:
    """
    Get upper and lower bound probability of detection given test/observations of each bait-prey combination
    """
    params:
        pseudo_n=config["pseudo_n"],
        id_pattern= config["id_pattern"]
    input:
        method_aggregate=f"work_folder{pn}/inferred_search_space/aggregated/methods/{{data}}_experimental_wise.csv"
    output:
        full_detection=f"work_folder{pn}/analysis/POD/{{network_type}}/POD_{{data}}.pq"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[f"{params.id_pattern}_bait"] != inferred_negative_df[f"{params.id_pattern}_prey"]
            ].copy()

        if wildcards.network_type == "undirectional":
            inferred_negative_df["id_var"] = inferred_negative_df[[f"{params.id_pattern}_bait", f"{params.id_pattern}_prey"]].apply(
                lambda x: "_".join(sorted(x)), axis=1)

            s = datetime.datetime.now()
            inferred_negative_df.sort_values("id_var", inplace=True)
            inferred_negative_mat = inferred_negative_df.to_numpy()
            aggregated_negative_mat = np.zeros_like(inferred_negative_mat)
            prev_bait, prev_prey, prev_n_observed, prev_n_tested, prev_pids, _, prev_id = inferred_negative_mat[0]
            for i in range(1, inferred_negative_mat.shape[0]):
                c_bait, c_prey, c_n_observed, c_n_tested, c_pid, _, c_id = inferred_negative_mat[i] # CLID is not used yet, but should be fixed later
                pids = set(c_pid.split(";"),)
                c_bait, c_prey = order_prot = sorted([c_bait, c_prey])
                
                if prev_id == c_id:
                    prev_n_observed += c_n_observed
                    prev_n_tested += c_n_tested
                    prev_pids |= pids 
                else:
                    aggregated_negative_mat[i-1] = [
                        prev_bait,
                        prev_prey,
                        prev_n_observed,
                        prev_n_tested,
                        ";".join(prev_pids),
                        "",
                        prev_id
                    ]
                    prev_bait, prev_prey, prev_n_observed, prev_n_tested, prev_pids, prev_id = c_bait, c_prey, c_n_observed, c_n_tested, pids, c_id
                
            aggregated_negative_mat[i] = [
                prev_bait,
                prev_prey,
                prev_n_observed,
                prev_n_tested,
                ";".join(prev_pids),
                "",
                prev_id
            ]

            ppis_joined_idx = aggregated_negative_mat[:, 1] != 0
            print(f"joined {-sum(ppis_joined_idx-1)} out of {inferred_negative_mat.shape[0]} rows in {(datetime.datetime.now() - s).total_seconds()} seconds", flush=True)
            aggregated_negative_mat = aggregated_negative_mat[ppis_joined_idx, :]
            undirectional_negative_df = pd.DataFrame(
                aggregated_negative_mat,
                columns=inferred_negative_df.columns
            )

            flipped_df = undirectional_negative_df.copy()
            flipped_df[["gene_name_bait", "gene_name_prey"]] = (
                flipped_df[["gene_name_prey", "gene_name_bait"]]
            )

            inferred_negative_df = (
                pd.concat([undirectional_negative_df, flipped_df], ignore_index=True)
            )

        inferred_negative_df["n_observed"] = inferred_negative_df["n_observed"].astype(int)
        inferred_negative_df["n_tested"] = inferred_negative_df["n_tested"].astype(int)

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
        inferred_negative_df["pair_id"] = range(inferred_negative_df.shape[0])

        inferred_negative_df.to_parquet(
            output.full_detection,
            index=False
        )


rule all_methods_filter_out_cell_line:
    """
    Get upper and lower bound probability of detection given test/observations of each bait-prey combination per cell line
    """
    params:
        pseudo_n=config["pseudo_n"],
        id_pattern= config["id_pattern"]
    input:
        method_aggregate=f"work_folder{pn}/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    output:
        full_detection=f"work_folder{pn}/analysis/POD/{{network_type}}/POD_cell_line.pq"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[f"{params.id_pattern}_bait"] != inferred_negative_df[f"{params.id_pattern}_prey"]
            ].copy()

        if wildcards.network_type == "undirectional":
            inferred_negative_df["id_var"] = inferred_negative_df[[f"{params.id_pattern}_bait", f"{params.id_pattern}_prey", "CVCL"]].apply(
                lambda x: "_".join(sorted(x[:2]))+"_" + x[3]  , axis=1)

            s = datetime.datetime.now()
            inferred_negative_df.sort_values("id_var", inplace=True)
            inferred_negative_mat = inferred_negative_df.to_numpy()
            aggregated_negative_mat = np.zeros_like(inferred_negative_mat)
            prev_bait, prev_prey, prev_n_observed, prev_n_tested, prev_pids, prev_cl, prev_id = inferred_negative_mat[0]
            for i in range(1, inferred_negative_mat.shape[0]):
                c_bait, c_prey, c_n_observed, c_n_tested, c_pid, c_cl, c_id = inferred_negative_mat[i]
                pids = set(c_pid.split(";"),)
                c_bait, c_prey = order_prot = sorted([c_bait, c_prey])
                
                if prev_id == c_id:
                    prev_n_observed += c_n_observed
                    prev_n_tested += c_n_tested
                    prev_pids |= pids 
                else:
                    aggregated_negative_mat[i-1] = [
                        prev_bait,
                        prev_prey,
                        prev_n_observed,
                        prev_n_tested,
                        ";".join(prev_pids),
                        prev_cl,
                        prev_id
                    ]
                    prev_bait, prev_prey, prev_n_observed, prev_n_tested, prev_pids, prev_cl, prev_id = c_bait, c_prey, c_n_observed, c_n_tested, pids, c_cl, c_id
                
            aggregated_negative_mat[i] = [
                prev_bait,
                prev_prey,
                prev_n_observed,
                prev_n_tested,
                ";".join(prev_pids),
                prev_cl,
                prev_id
            ]

            ppis_joined_idx = aggregated_negative_mat[:, 1] != 0
            print(f"joined {-sum(ppis_joined_idx-1)} out of {inferred_negative_mat.shape[0]} rows in {(datetime.datetime.now() - s).total_seconds()} seconds", flush=True)
            aggregated_negative_mat = aggregated_negative_mat[ppis_joined_idx, :]
            undirectional_negative_df = pd.DataFrame(
                aggregated_negative_mat,
                columns=inferred_negative_df.columns
            )

            flipped_df = undirectional_negative_df.copy()
            flipped_df[["gene_name_bait", "gene_name_prey"]] = (
                flipped_df[["gene_name_prey", "gene_name_bait"]]
            )

            inferred_negative_df = (
                pd.concat([undirectional_negative_df, flipped_df], ignore_index=True)
            )

        inferred_negative_df["n_observed"] = inferred_negative_df["n_observed"].astype(int)
        inferred_negative_df["n_tested"] = inferred_negative_df["n_tested"].astype(int)

        global_pod_by_cvcl = (
            inferred_negative_df
                .groupby("CVCL")
                .apply(lambda x: x["n_observed"].sum() / x["n_tested"].sum
                )
                )
        inferred_negative_df["alpha_prior"] = inferred_negative_df["CVCL"].map(global_pod_by_cvcl)
        inferred_negative_df["alpha_prior"] = inferred_negative_df["alpha_prior"] * params.pseudo_n
        inferred_negative_df["beta_prior"] = params.pseudo_n - inferred_negative_df["alpha_prior"]
        
        inferred_negative_df["alpha_post"] = prior_alpha + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"] = prior_beta + inferred_negative_df["n_tested"] - inferred_negative_df["n_observed"]

        inferred_negative_df["p"] = inferred_negative_df["alpha_post"] / (
                inferred_negative_df["alpha_post"] + inferred_negative_df["beta_post"])

        inferred_negative_df["lower_bound_pod"] = beta.ppf(0.025,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["upper_bound_pod"] = beta.ppf(0.975,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["pair_id"] = range(inferred_negative_df.shape[0])

        inferred_negative_df.to_parquet(
            output.full_detection,
            index=False
        )



def t_threshold_degree(df, t, greater=True, n = 5):
    if greater:
        df_t = df[df["lower_bound_pod"] > t]
    else:
        df_t = df[(df["n_observed"] == 0) & (df["n_tested"] >= n)]

    df_bait_degree = df_t.groupby("gene_name_bait",as_index=False).size()
    df_bait_degree = df_bait_degree.rename({
        "gene_name_bait": "gene_name",
        "size": "degree_bait"
    },axis=1)
    df_prey_degree = df_t.groupby("gene_name_prey",as_index=False).size()
    df_prey_degree = df_prey_degree.rename({
        "gene_name_prey": "gene_name",
        "size": "degree_prey"
    },axis=1)
    t_degree = df_bait_degree.merge(df_prey_degree,on="gene_name",how="outer").fillna(0)

    return t_degree

rule get_hcl_degree:
    input:
        ppis = f"work_folder{pn}/analysis/POD/POD_{{method}}.csv"
    output:
        pos_01 = f"work_folder{pn}/degree/{{method}}_t0.1.csv",
        neg_gt5 = f"work_folder{pn}/degree/{{method}}_gt_5.csv"
    run:
        df = pd.read_csv(input.ppis, sep="\t")
        t_threshold_degree(df, 0.1).to_csv(output.pos_01, sep="\t", index=False)
        t_threshold_degree(df,0.1, greater=False).to_csv(output.neg_gt5, sep="\t", index=False)