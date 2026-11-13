from scipy.stats import spearmanr

rule get_all_balance_metrics:
    input:
        metrics=expand(
            "work_folder/subsets/train/permuted/{permutation}/balance/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_degree.csv",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["classification_sources"]["pod"]["positive"]["limits"],
            neg_limit=config["classification_sources"]["pod"]["negative"]["limits"],
            random=["", "-random"],
        ),
    output:
        all_models="work_folder/analysis/subset_similarity/balance/all_metrics_{network_type}.csv",
    shell:
        """
        echo "permutation\tdataset\tpos_limit\tneg_limit\trandom\tbait_degree_delta\tprey_degree_delta\tspearman_bait\tspearman_prey\tnum_edges" > {output.all_models}
        cat {input.metrics} >> {output.all_models}
        """


rule plot_degree_balance:
    input:
        metrics="work_folder/analysis/subset_similarity/balance/all_metrics_{network_type}.csv",
    output:
        plot="work_folder/analysis/subset_similarity/balance/plots/{network_type}_degree_balance.png",
    log:
        "logs/analysis/subset_similarity/balance/plots/{network_type}_degree_balance.log",
    script:
        "scripts/plot_degree_balance.py"



rule get_train_degree_balance:
    input:
        train_pos="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
    output:
        balance="work_folder/subsets/train/permuted/{permutation}/balance/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_degree.csv",
    run:
        pos_df = pd.read_csv(input.train_pos, sep="\t")
        neg_df = pd.read_csv(input.train_neg, sep="\t")

        if wildcards.network_type == "directional":
            graph_type = nx.DiGraph
        elif wildcards.network_type == "undirectional":
            graph_type = nx.Graph
        else:
            raise ValueError(f"{wildcards.network_type} is not a valid network type.")

        pos_G = nx.from_pandas_edgelist(pos_df, "bait", "prey", create_using=graph_type)
        neg_G = nx.from_pandas_edgelist(neg_df, "bait", "prey", create_using=graph_type)
        n_edges = pos_G.number_of_edges() + neg_G.number_of_edges()

        if wildcards.network_type == "directional":
            pos_out = dict(pos_G.out_degree())
            neg_out = dict(neg_G.out_degree())
            pos_in = dict(pos_G.in_degree())
            neg_in = dict(neg_G.in_degree())

            all_bait = set(pos_out) | set(neg_out)
            all_prey = set(pos_in) | set(neg_in)
            degree_bait_delta = sum(abs(pos_out.get(n, 0) - neg_out.get(n, 0)) for n in all_bait)
            degree_prey_delta = sum(abs(pos_in.get(n, 0) - neg_in.get(n, 0)) for n in all_prey)

            bait_pos_deg = [pos_out.get(n, 0) for n in all_bait]
            bait_neg_deg = [neg_out.get(n, 0) for n in all_bait]
            prey_pos_deg = [pos_in.get(n, 0) for n in all_prey]
            prey_neg_deg = [neg_in.get(n, 0) for n in all_prey]
            spearman_bait = spearmanr(bait_pos_deg, bait_neg_deg).correlation if len(all_bait) > 1 else float("nan")
            spearman_prey = spearmanr(prey_pos_deg, prey_neg_deg).correlation if len(all_prey) > 1 else float("nan")
        else:
            pos_deg = dict(pos_G.degree())
            neg_deg = dict(neg_G.degree())
            all_nodes = set(pos_deg) | set(neg_deg)

            degree_bait_delta = sum(abs(pos_deg.get(n, 0) - neg_deg.get(n, 0)) for n in all_nodes)
            degree_prey_delta = 0  # avoid double-counting the same delta in bait+prey sums downstream

            pos_vals = [pos_deg.get(n, 0) for n in all_nodes]
            neg_vals = [neg_deg.get(n, 0) for n in all_nodes]
            spearman_bait = spearmanr(pos_vals, neg_vals).correlation if len(all_nodes) > 1 else float("nan")
            spearman_prey = spearman_bait

        with open(output.balance, "w") as w:
            w.write(
                f"{wildcards.permutation}\t{wildcards.dataset}\t{wildcards.pos_limit}\t{wildcards.neg_limit}\t{wildcards.random!= ""}\t"
                f"{degree_bait_delta}\t{degree_prey_delta}\t{spearman_bait}\t{spearman_prey}\t{n_edges}\n"
            )

rule define_balanced_negative:
    input:
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
    output:
        selected_negative="work_folder/subsets/{source}/test/{dataset}_{network_type}_selected_hrni.tsv",
    run:
        df_pos = pd.read_csv(input.test_pos, sep="\t")[["bait", "prey"]]
        df_negative = pd.read_csv(input.test_neg, sep="\t")[["bait", "prey"]]
        if df_negative.shape[0] > df_pos.shape[0]:
            df_negative = df_negative.sample(df_pos.shape[0], random_state=1234)
        df_negative.to_csv(output.selected_negative, sep="\t", index=False)


rule check_against_all_non_observed:
    input:
        selected_negative="work_folder/subsets/{source}/test/{dataset}_{network_type}_selected_hrni.tsv",
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        saved_model="work_folder/classification/{source}/{classifier}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        no_selected_negative="work_folder/subsets/{source}/test/no/{dataset}_{network_type}_selected_no.tsv"
    script:
        "scripts/check_against_all_no.py"


rule get_model_metrics:
    input:
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
        selected_negative="work_folder/subsets/{source}/test/{dataset}_{network_type}_selected_{neg_type}.tsv",
        saved_model="work_folder/classification/{source}/{classifier}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        metrics="work_folder/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics_{neg_type}.txt",
        roc_png="work_folder/classification/{source}/{classifier}/metrics/plot/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_roc_curve_{neg_type}.png",
    log:
        "logs/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_metrics_{neg_type}.log",
    threads: 1
    resources:
        mem_gb=3,
    params:
        script_location="src/PPIClassification/ModelEvaluation/evaluate_model.py",
    shell:
        """
        python3 {params.script_location} \
            --pos_data_file {input.test_pos} \
            --neg_data_file {input.test_neg} \
            --protein_embeddings_file {input.protein_embeddings} \
            --model_file {input.saved_model} \
            --output_file {output.metrics} \
            --plot_roc_png {output.roc_png} \
            --network_type  {wildcards.network_type} \
            --neg_input_file {input.selected_negative}> {log} 2>&1
        """


def source_grid_metrics(input_path, wildcards):
    source = config["classification_sources"][wildcards.source]
    return expand(
        input_path,
        **dict(wildcards.items()),  # source, network_type, esm_model (+ classifier if the rule has it)
        dataset=source["datasets"],
        pos_limit=source["positive"]["limits"],
        neg_limit=source["negative"]["limits"],
        random=["", "-random"],
        negative=["hrni", "no"]
    )


rule model_all_metrics:
    input:
        metrics=lambda wc: source_grid_metrics(
            "work_folder/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics_{negative}.txt",
            wc
        ),
    output:
        all_models="work_folder/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.csv",
    log:
        "logs/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.log",
    run:
        with open(output[0], "w") as w:
            w.write("model\troc_auc\tacc_I\tacc_NI\tsamples\n")
            for metric_file in input.metrics:
                with open(metric_file, "r") as f:
                    line_out = [line.strip().split(": ")[1] for line in f]
                    line_out = "\t".join(line_out)
                    model_name = metric_file.split("/")[-1].replace(
                        "_metrics.txt", ""
                    )
                    w.write(
                        model_name + "\t" + line_out + "\n"
                    )


def all_metrics_in_run(input_path, wildcards):
    source = config["classification_sources"]
    classification_targets = config["classification_targets"]
    return expand(
        input_path,
        source = classification_targets["source"],
        classifier = classification_targets["classifier"],
        network_type=classification_targets["network_type"],
        esm_model=classification_targets["esm_model"]
    )


rule plot_hrni_no_auc_arrows:
    input:
        metrics=lambda wc: all_metrics_in_run(
            "work_folder/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.csv",
            wc
        ),
    output:
        plot="work_folder/classification/plots/hrni_no_auc_arrows.png",
    script:
        "scripts/plot_hrni_no_auc_arrows.R"


rule within_pair_similarity_auc:
    input:
        train_pos="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_hrni="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_neg.csv",
        train_no="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-random_neg.csv",
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_hrni="work_folder/subsets/{source}/test/{dataset}_{network_type}_selected_hrni.tsv",
        test_no="work_folder/subsets/{source}/test/{dataset}_{network_type}_selected_no.tsv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        auc="work_folder/classification/{source}/embedding_sim_auc/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_{esm_model}.tsv",
    resources:
        mem_gb=8,
    script:
        "scripts/within_pair_similarity_auc.py"
        



rule collect_all_cosine_similarites:
    input:
        sims = lambda wc: source_grid_metrics(
            "work_folder/classification/{source}/embedding_sim_auc/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_{esm_model}.tsv",
            wc
        )
    output:
        auc="work_folder/classification/{source}/embedding_sim_auc/all_{network_type}_{esm_model}.tsv",
    shell:
        "awk 'FNR > 1 || NR == 1' {input.sims} > {output.auc}"


rule collect_all_similarities:
    input:
        all_sims = lambda wc: all_metrics_in_run(
            "work_folder/classification/{source}/embedding_sim_auc/all_{network_type}_{esm_model}.tsv",
            wc),
    output:
        all_auc = "work_folder/classification/embedding_sim_auc/all_classifications.tsv"
    shell:
        "awk 'FNR > 1 || NR == 1' {input.all_sims} > {output.all_auc}"


rule hrni_no_delta_auc:
    input:
        metrics="work_folder/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.csv",
    output:
        by_config="work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_{esm_model}_delta_auc_by_config.tsv",
        summary="work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_{esm_model}_delta_auc_summary.tsv",
    log:
        "logs/classification/{source}/{classifier}/hrni_no_delta/{network_type}_{esm_model}_delta_auc.log",
    run:
        import re
        from scipy.stats import ttest_1samp

        model_re = re.compile(
            r"^(?P<dataset>[a-zA-Z0-9]+)_[a-zA-Z]+_limit_(?P<neg_limit>[0-9.]+)"
            r"_poslim_(?P<pos_limit>[0-9.]+|all)(?P<random>-random)?_model_[A-Za-z0-9]+"
            r"_metrics_(?P<test_neg>hrni|no)\.txt$"
        )
        df = pd.read_csv(input.metrics, sep="\t")
        df = pd.concat([df, df["model"].str.extract(model_re)], axis=1)
        assert df["dataset"].notna().all(), f"unparsed model names: {df.loc[df['dataset'].isna(), 'model'].tolist()}"
        df = df[df["test_neg"] == "hrni"]  # always score against the HRNI test set
        df["train_type"] = df["random"].notna().map({True: "no", False: "hrni"})

        pivot = df.pivot_table(
            index=["dataset", "pos_limit", "neg_limit"],
            columns="train_type", values="roc_auc",
        )
        pivot["delta_auc"] = pivot["hrni"] - pivot["no"]
        pivot = pivot.reset_index()
        pivot.to_csv(output.by_config, sep="\t", index=False)

        rows = []
        for dataset, sub in pivot.groupby("dataset"):
            t, p = ttest_1samp(sub["delta_auc"], 0.0) if len(sub) > 1 else (float("nan"), float("nan"))
            rows.append({
                "dataset": dataset,
                "n_configs": len(sub),  # n(pos_limit x neg_limit)
                "mean_delta_auc": sub["delta_auc"].mean(),
                "sem": sub["delta_auc"].std(ddof=1) / len(sub) ** 0.5,
                "t": t,
                "p": p,
            })
        pd.DataFrame(rows).to_csv(output.summary, sep="\t", index=False)


rule plot_hrni_no_delta_auc:
    input:
        by_config="work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_{esm_model}_delta_auc_by_config.tsv",
    output:
        plot="work_folder/classification/{source}/{classifier}/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.png",
    log:
        "logs/classification/{source}/{classifier}/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.log",
    script:
        "scripts/plot_hrni_no_delta_auc.py"


rule plot_hrni_no_delta_auc_prost_esm2:
    input:
        esm_by_config="work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_ESM2_delta_auc_by_config.tsv",
        prost_by_config="work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_ProstT5_delta_auc_by_config.tsv",
    output:
        plot="work_folder/classification/{source}/{classifier}/hrni_no_delta/plots/{network_type}_embedding_delta_auc_boxplot.png",
    log:
        "logs/classification/{source}/{classifier}/hrni_no_delta/plots/{network_type}_embeddings_delta_auc_boxplot.log",
    script:
        "scripts/plot_hrni_no_delta_auc_prost.py"

