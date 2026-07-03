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
        method_aggregate="work_folder/inferred_search_space/aggregated/methods/{data}_experimental_wise.csv",
        fasta="work_folder/protein_sequences/uniprot_canonical.fasta"
    output:
        full_detection="work_folder/analysis/POD/{network_type}/POD_{data}.pq"
    log:
        "logs/analysis/POD/{network_type}/POD_{data}.log"
    run:
        bait_col = f"{params.id_pattern}_bait"
        prey_col = f"{params.id_pattern}_prey"

        fasta_ids = set()
        with open(input.fasta) as _fasta:
            for line in _fasta:
                if line.startswith(">"):
                    parts = line.split("|")
                    if len(parts) >= 2:
                        fasta_ids.add(parts[1])

        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t",
            dtype={
                bait_col: "string",
                prey_col: "string"
            }
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[bait_col].isin(fasta_ids) &
            inferred_negative_df[prey_col].isin(fasta_ids)
        ].copy()

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[bait_col] != inferred_negative_df[prey_col]
            ].copy()

        if wildcards.network_type == "undirectional":
            inferred_negative_df["id_var"] = inferred_negative_df[[bait_col, prey_col]].apply(
                lambda x: "_".join(sorted(x)), axis=1)

            s = datetime.datetime.now()
            n_before = len(inferred_negative_df)
            inferred_negative_df[[bait_col, prey_col]] = pd.DataFrame(
                inferred_negative_df[[bait_col, prey_col]].apply(sorted, axis=1).tolist(),
                index=inferred_negative_df.index,
                columns=[bait_col, prey_col]
            )
            orig_cols = inferred_negative_df.columns.tolist()
            inferred_negative_df = (
                inferred_negative_df
                .groupby("id_var", sort=False)
                .agg({
                    bait_col: "first",
                    prey_col: "first",
                    "n_tested": "sum",
                    "n_observed": "sum",
                    "pubmed_id": lambda x: ";".join({p for ids in x for p in ids.split(";")}),
                    "CVCL": "first",
                })
                .reset_index()
            )[orig_cols]
            print(f"joined {n_before - len(inferred_negative_df)} out of {n_before} rows in {(datetime.datetime.now() - s).total_seconds()} seconds", flush=True)

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
        method_aggregate="work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    output:
        full_detection="work_folder/analysis/POD/{network_type}/POD_cell_line.pq"
    log:
        "logs/analysis/POD/{network_type}/POD_cell_line.log"
    run:
        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t"
        )

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[bait_col] != inferred_negative_df[prey_col]
            ].copy()

        if wildcards.network_type == "undirectional":
            prot_a = inferred_negative_df[[bait_col, prey_col]].min(axis=1)
            prot_b = inferred_negative_df[[bait_col, prey_col]].max(axis=1)
            inferred_negative_df["id_var"] = prot_a + "_" + prot_b + "_" + inferred_negative_df["CVCL"]
            
            inferred_negative_df.sort_values("id_var", inplace=True)
            inferred_negative_mat = inferred_negative_df.to_numpy()
            aggregated_negative_mat = np.zeros_like(inferred_negative_mat)
            prev_bait, prev_prey, prev_n_observed, prev_n_tested, prev_pids, prev_cl, prev_id = inferred_negative_mat[0]
            prev_pids = set(prev_pids.split(";"))
            s = datetime.datetime.now()
            for i in range(1, inferred_negative_mat.shape[0]):
                c_bait, c_prey, c_n_observed, c_n_tested, c_pid, c_cl, c_id = inferred_negative_mat[i]
                pids = set(c_pid.split(";"))
                c_bait, c_prey = sorted([c_bait, c_prey])
                
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
            inferred_negative_df = pd.DataFrame(
                aggregated_negative_mat,
                columns=inferred_negative_df.columns
            )

        inferred_negative_df["n_observed"] = inferred_negative_df["n_observed"].astype(int)
        inferred_negative_df["n_tested"] = inferred_negative_df["n_tested"].astype(int)

        pod_by_cvcl = (
            inferred_negative_df
                .groupby("CVCL")[["n_observed", "n_tested"]]
                .sum()
        )
        pod_by_cvcl["cl_pod"] = (
            pod_by_cvcl["n_observed"]/pod_by_cvcl["n_tested"]
        )

        inferred_negative_df["alpha_prior"] = inferred_negative_df["CVCL"].map(pod_by_cvcl["cl_pod"])
        inferred_negative_df["alpha_prior"] = inferred_negative_df["alpha_prior"] * params.pseudo_n
        inferred_negative_df["beta_prior"] = params.pseudo_n - inferred_negative_df["alpha_prior"]
        
        inferred_negative_df["alpha_post"] = inferred_negative_df["alpha_prior"] + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"] = inferred_negative_df["beta_prior"] + inferred_negative_df["n_tested"] - inferred_negative_df["n_observed"]

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

