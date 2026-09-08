rule get_other_ppis:
    input:
        miTab = "work_folder/data/intact/human.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        other_ppis = "work_folder/formated/other_method_ppis.csv"
    log:
        "logs/formated/other_method_ppis.log"
    run:
        mitab_df = filter_mitab(input.miTab)
        mitab_df["detection_method"] = mitab_df["detection_method"].str.replace(":", "-")

        selected_methods = set(config["ms"]) | set(config["y2h"])
        other_df = mitab_df[~mitab_df["detection_method"].isin(selected_methods)]

        other_df = other_df[["IDA", "IDB", "detection_method"]].dropna()
        other_df = other_df.rename(columns={"IDA": "prot_a", "IDB": "prot_b"})

        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        other_df = other_df.merge(gene_name_df, left_on="prot_a", right_on="uniprot_id")
        del other_df["uniprot_id"]
        other_df = other_df.merge(gene_name_df, left_on="prot_b", right_on="uniprot_id", suffixes=("_a", "_b"))
        del other_df["uniprot_id"]

        other_df = other_df.drop_duplicates()
        other_df = other_df[other_df["prot_a"] != other_df["prot_b"]] # drop homomers
        other_df.to_csv(output.other_ppis, sep="\t", index=None)


rule contradiciton_rate:
    params:
        HRNI_limits = [1,3,5],
        interaction_limit = 1
    input:
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq",
        other_ppi_intact = "work_folder/formated/other_method_ppis.csv",
        method_groupings = "data/method_grouping.yaml"
    output:
        detection_statistics = "work_folder/analysis/other_methods/detection_stats/{dataset}_stats.csv"
    script:
        "scripts/method_agreement.py"


rule plot_contradiction_rate:
    input:
        stats_file = "work_folder/analysis/other_methods/detection_stats/{dataset}_stats.csv"
    output:
        png = "work_folder/analysis/other_methods/detection_stats/plot/{dataset}_contradiction_rate.png"
    script:
        "scripts/plot_contradiction_rate.py"


rule get_fdr_for:
    params:
        hrni_limit = [1,3],
        hri_limit = 0.15
    input:
        ms_pod = "work_folder/analysis/POD/undirectional/POD_ms.pq",
        y2h_pod = "work_folder/analysis/POD/undirectional/POD_y2h.pq",
        y2h_leave_out = "work_folder/inferred_search_space/experimental_method/20211142_MI-0397.csv",
        ms_leave_out = "work_folder/inferred_search_space/experimental_method/32707033_MI-0096.csv"
    output:
        rates = "work_folder/analysis/other_methods/leave_out/fdr_for.csv"
    log:
        "logs/analysis/other_methods/leave_out_fdr_for.log"
    threads:
        1
    run:
        import pyarrow.parquet as pq
        from scipy.stats import beta

        id_pattern = config["id_pattern"]
        bait_col, prey_col = f"{id_pattern}_bait", f"{id_pattern}_prey"

        def evaluate(pod_path, leave_out_path, dataset):
            # 1. load pod (arrow reader, the POD tables run to ~1e8 rows)
            pod = pq.read_table(pod_path, columns=[
                bait_col, prey_col, "n_tested", "n_observed",
                "alpha_post", "beta_post"],
            ).to_pandas(split_blocks=True, self_destruct=True)

            pod["pair_id"] = pod[[bait_col, prey_col]].apply(lambda row: ":".join(sorted(row)), axis=1)
            
            row = pod.iloc[0]
            prior_alpha = row["alpha_post"] - row["n_observed"]
            prior_beta = row["beta_post"] - (row["n_tested"] - row["n_observed"])

            # 2. remove protein pairs from pod from correct leave out.
            held = pd.read_csv(leave_out_path, sep="\t")
            held["pair_id"] = held[[bait_col, prey_col]].apply(lambda row: ":".join(sorted(row)), axis=1)
            held = held.groupby("pair_id", as_index=False).agg(
                held_tested=("n_tested", "sum"), held_observed=("n_observed", "sum"))

            df = pod.merge(held, on="pair_id", how="inner")
            df["n_tested_loo"] = (df["n_tested"] - df["held_tested"]).clip(lower=0)
            df["n_observed_loo"] = (df["n_observed"] - df["held_observed"]).clip(lower=0)            

            # 3. calculate the Q.025 for affected rows
            df["alpha_post_loo"] = prior_alpha + df["n_observed_loo"]
            df["beta_post_loo"] = prior_beta + df["n_tested_loo"] - df["n_observed_loo"]
            df["lower_bound_pod_loo"] = beta.ppf(
                0.025, df["alpha_post_loo"], df["beta_post_loo"])
            
            #df["hri"] = df["lower_bound_pod_loo"] > params.hri_limit
            df["hri"] = df["n_observed"] != 0

            fp_dict = []
            for hrni_lim in params.hrni_limit:
                df["hrni"] = (df["n_tested_loo"] >= hrni_lim) & (df["n_observed_loo"] == 0) 
                ss_df = df[
                    (df["hri"]) |
                    (df["hrni"])
                ]
                # unit is one held-out test; "positive" means that test detected the pair
                hri_rows, hrni_rows = ss_df[ss_df["hri"]], ss_df[ss_df["hrni"]]
                tp = hri_rows["held_observed"].sum()
                fn = hri_rows["held_tested"].sum() - tp
                fp = hrni_rows["held_observed"].sum()
                tn = hrni_rows["held_tested"].sum() - fn
                FDR = fp/(fp+tp)
                FOR = fn/(tn+fn)
                cov = ss_df.shape[0]
                fp_dict.append([hrni_lim, tp, fp, fn, tn, FDR, FOR, cov])

            return fp_dict

        results = []
        for dataset, pod_path, held_path in [
            ("y2h", input.y2h_pod, input.y2h_leave_out),
            ("ms", input.ms_pod, input.ms_leave_out),
        ]:
            for row in evaluate(pod_path, held_path, dataset):
                print(dataset, row, flush=True)
                results.append([dataset] + row)

        pd.DataFrame(results, columns=[
            "dataset", "hrni_limit", "TP", "FP", "FN", "TN",
            "fdr", "for", "n_called",
        ]).to_csv(output.rates, sep="\t", index=False)
