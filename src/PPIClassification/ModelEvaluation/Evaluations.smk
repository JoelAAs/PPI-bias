from scipy.stats import spearmanr

rule plot_degree_balance:
    input:
        metrics="work_folder/analysis/subset_similarity/balance/all_metrics_{network_type}.csv",
    output:
        plot="work_folder/analysis/subset_similarity/balance/plots/{network_type}_degree_balance.png",
    log:
        "logs/analysis/subset_similarity/balance/plots/{network_type}_degree_balance.log",
    script:
        "scripts/plot_degree_balance.py"


rule get_non_observbed_negatives:
    input:
        test_neg="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
    output:
        no_selected_negative="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg_no.csv"
    script:
        "scripts/get_equivalent_negative.py"


rule get_model_metrics:
    input:
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
        test_neg_no="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg_no.csv",
        saved_model="work_folder/classification/{source}/{classifier}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
        train_pos="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        metrics="work_folder/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics.txt",
        roc_png="work_folder/classification/{source}/{classifier}/metrics/plot/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_roc_curve.png",
        roc_png_no="work_folder/classification/{source}/{classifier}/metrics/plot/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_roc_curve_no.png",
        pred="work_folder/classification/{source}/{classifier}/predictions/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_pred.tsv.gz",
        pred_no="work_folder/classification/{source}/{classifier}/predictions/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_pred_no.tsv.gz",
    log:
        "logs/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_metrics.log",
    threads: 1
    resources:
        mem_gb=3,
    script:
        "evaluate_model.py"


def source_grid_metrics(input_path, wildcards):
    source = config["classification_sources"][wildcards.source]
    return expand(
        input_path,
        **dict(wildcards.items()),  # source, network_type, esm_model (+ classifier if the rule has it)
        dataset=source["datasets"],
        pos_limit=source["positive"]["limits"],
        neg_limit=source["negative"]["limits"],
        random=["", "-random"]
    )


rule model_all_metrics:
    input:
        metrics=lambda wc: source_grid_metrics(
            "work_folder/classification/{source}/{classifier}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics.txt",
            wc
        ),
    output:
        all_models="work_folder/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.csv",
    log:
        "logs/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.log",
    run:
        with open(output[0], "w") as w:
            w.write("model\troc_auc\tacc_I\tacc_NI\tsamples\troc_auc_no\tacc_NI_no\ttrain_samples\tsamples_no\troc_auc_adj\troc_auc_adj_no\n")
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

rule full_run_metrics:
    input:
        metrics=lambda wc: all_metrics_in_run(
            "work_folder/classification/{source}/{classifier}/all_metrics_{network_type}_{esm_model}.csv",
            wc
        ),
    output:
        metrics = "work_folder/classification/full_run_metrics.tsv"
    run:
        dfs = []
        for path in input.metrics:
            source, classifier = path.split("/")[-3:-1]
            df = pd.read_csv(path, sep="\t")
            df.insert(0, "classifier", classifier)
            df.insert(0, "source", source)
            dfs.append(df)
        pd.concat(dfs).to_csv(output.metrics, sep="\t", index=False)

rule plot_hrni_no_auc_arrows:
    input:
        metrics = "work_folder/classification/full_run_metrics.tsv"
    output:
        plot="work_folder/classification/plots/hrni_no_auc_arrows.png",
    script:
        "scripts/plot_hrni_no_auc_arrows.R"





# Pair-wise Embedding similarity
rule within_pair_similarity_auc:
    input:
        train_pos="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_hrni="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_neg.csv",
        train_no="work_folder/subsets/{source}/train/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-random_neg.csv",
        test_pos="work_folder/subsets/{source}/test/{dataset}_{network_type}_pos.csv",
        test_hrni="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg.csv",
        test_no="work_folder/subsets/{source}/test/{dataset}_{network_type}_neg_no.csv",
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


## Similarity adjusted prediction

# rule get_pair_wise_similarity:
#     input:
#         metrics = "work_folder/classification/full_run_metrics.tsv"
#         sim_auc = "work_folder/classification/embedding_sim_auc/all_classifications.tsv"
#     output:
#         delta_metrics = "work_folder/classification/full_run_metrics_sim.tsv"
#     run:



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
            r"_poslim_(?P<pos_limit>[0-9.]+|all)(?P<random>-random)?_model_[A-Za-z0-9]+$"
        )
        df = pd.read_csv(input.metrics, sep="\t")
        df = pd.concat([df, df["model"].str.extract(model_re)], axis=1)
        assert df["dataset"].notna().all(), f"unparsed model names: {df.loc[df['dataset'].isna(), 'model'].tolist()}"
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

