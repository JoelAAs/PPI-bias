from scipy.stats import spearmanr


rule get_all_balance_metrics:
    input:
        metrics=expand(
            "work_folder/subsets/train/permuted/{permutation}/balance/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_degree.csv",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
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
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/test/{dataset}_{network_type}_neg.csv",
    output:
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_hrni.tsv",
    run:
        df_pos = pd.read_csv(input.test_pos, sep="\t")[["bait", "prey"]]
        df_negative = pd.read_csv(input.test_neg, sep="\t")[["bait", "prey"]]
        if df_negative.shape[0] > df_pos.shape[0]:
            df_negative = df_negative.sample(df_pos.shape[0], random_state=1234)
        df_negative.to_csv(output.selected_negative, sep="\t", index=False)


rule define_no_negative_test_set:
    # Fixed per (dataset, network_type), same as define_balanced_negative above: used as the
    # shared "NO" (non-observed) negative test set by every classifier/permutation/pos_limit/
    # neg_limit/esm_model model, so it's computed once rather than once per model.
    input:
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_hrni.tsv",
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv"
    output:
        no_selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_no.tsv"
    script:
        "scripts/get_equivalent_negative.py"


rule get_model_metrics_permuted_hnri:
    input:
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/test/{dataset}_{network_type}_neg.csv",
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_hrni.tsv",
        saved_model="work_folder/classification/{classifier}/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        metrics="work_folder/classification/{classifier}/permuted/{permutation}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics_hrni.txt",
        roc_png="work_folder/classification/{classifier}/permuted/{permutation}/metrics/plot/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_roc_curve_hrni.png",
    log:
        "logs/classification/{classifier}/permuted/{permutation}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_metrics_hrni.log",
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


rule get_model_metrics_permuted_no:
    input:
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv",
        test_neg="work_folder/subsets/test/{dataset}_{network_type}_neg.csv",
        saved_model="work_folder/classification/{classifier}/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_no.tsv",
    output:
        metrics="work_folder/classification/{classifier}/permuted/{permutation}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_metrics_no.txt",
        roc_png="work_folder/classification/{classifier}/permuted/{permutation}/metrics/plot/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_roc_curve_no.png",
    log:
        "logs/classification/{classifier}/permuted/{permutation}/metrics/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_metrics_no.log",
    threads: 1
    resources:
        mem_gb=2,
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


rule all_metrics_permuted:
    input:
        metrics=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{{esm_model}}_metrics_{negative}.txt",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
            random=["", "-random"],
            negative=["hrni", "no"]
        ),
    output:
        all_models="work_folder/classification/{classifier}/permuted/all_metrics_{network_type}_{esm_model}.csv",
    log:
        "logs/classification/{classifier}/permuted/all_metrics_{network_type}_{esm_model}.log",
    run:
        with open(output[0], "w") as w:
            w.write("permutation\tmodel\troc_auc\tacc_I\tacc_NI\tsamples\n")
            for metric_file in input.metrics:
                with open(metric_file, "r") as f:
                    line_out = [line.strip().split(": ")[1] for line in f]
                    line_out = "\t".join(line_out)
                    model_name = metric_file.split("/")[-1].replace(
                        "_metrics.txt", ""
                    )
                    permutation = metric_file.split("/permuted/")[1].split("/")[0]
                    w.write(
                        permutation + "\t" + model_name + "\t" + line_out + "\n"
                    )


rule hrni_no_delta_auc:
    input:
        metrics="work_folder/classification/{classifier}/permuted/all_metrics_{network_type}_{esm_model}.csv",
    output:
        by_config="work_folder/classification/{classifier}/permuted/hrni_no_delta/{network_type}_{esm_model}_delta_auc_by_config.tsv",
        summary="work_folder/classification/{classifier}/permuted/hrni_no_delta/{network_type}_{esm_model}_delta_auc_summary.tsv",
    log:
        "logs/classification/{classifier}/permuted/hrni_no_delta/{network_type}_{esm_model}_delta_auc.log",
    run:
        import re
        from scipy.stats import ttest_1samp

        model_re = re.compile(
            r"^(?P<dataset>[a-zA-Z0-9]+)_[a-zA-Z]+_limit_(?P<neg_limit>[0-9.]+)"
            r"_poslim_(?P<pos_limit>[0-9.]+|all)(?P<random>-random)?_model_[A-Z0-9]+"
            r"_metrics_(?P<test_neg>hrni|no)\.txt$"
        )
        df = pd.read_csv(input.metrics, sep="\t")
        df = pd.concat([df, df["model"].str.extract(model_re)], axis=1)
        df = df[df["test_neg"] == "hrni"]  # always score against the HRNI test set
        df["train_type"] = df["random"].notna().map({True: "no", False: "hrni"})

        pivot = df.pivot_table(
            index=["dataset", "pos_limit", "neg_limit", "permutation"],
            columns="train_type", values="roc_auc",
        )
        pivot["delta_auc"] = pivot["hrni"] - pivot["no"]
        pivot = pivot.reset_index()
        pivot.to_csv(output.by_config, sep="\t", index=False)

        rows = []
        for dataset, sub in pivot.groupby("dataset"):
            t, p = ttest_1samp(sub["delta_auc"], 0.0)
            rows.append({
                "dataset": dataset,
                "n_configs": len(sub),  # n_permutations * n(pos_limit x neg_limit)
                "mean_delta_auc": sub["delta_auc"].mean(),
                "sem": sub["delta_auc"].std(ddof=1) / len(sub) ** 0.5,
                "t": t,
                "p": p,
            })
        pd.DataFrame(rows).to_csv(output.summary, sep="\t", index=False)


rule plot_hrni_no_delta_auc:
    input:
        by_config="work_folder/classification/{classifier}/permuted/hrni_no_delta/{network_type}_{esm_model}_delta_auc_by_config.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.png",
    log:
        "logs/classification/{classifier}/permuted/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.log",
    script:
        "scripts/plot_hrni_no_delta_auc.py"


rule pos_neg_limit_auc:
    input:
        metrics="work_folder/classification/{classifier}/permuted/all_metrics_{network_type}_{esm_model}.csv",
    output:
        by_config="work_folder/classification/{classifier}/permuted/pos_neg_limit_auc/{network_type}_{esm_model}_auc_by_config.tsv",
    log:
        "logs/classification/{classifier}/permuted/pos_neg_limit_auc/{network_type}_{esm_model}_auc_by_config.log",
    run:
        import re

        model_re = re.compile(
            r"^(?P<dataset>[a-zA-Z0-9]+)_[a-zA-Z]+_limit_(?P<neg_limit>[0-9.]+)"
            r"_poslim_(?P<pos_limit>[0-9.]+|all)(?P<random>-random)?_model_[A-Z0-9]+"
            r"_metrics_(?P<test_neg>hrni|no)\.txt$"
        )
        df = pd.read_csv(input.metrics, sep="\t")
        df = pd.concat([df, df["model"].str.extract(model_re)], axis=1)
        df = df[df["test_neg"] == "hrni"]  # always score against the HRNI test set
        df = df[df["random"].isna()]  # HRNI-trained models only, drop the "-random" (NO) models

        by_config = df.groupby(["dataset", "pos_limit", "neg_limit"])["roc_auc"].agg(
            mean_auc="mean", std_auc="std", n_permutations="count"
        ).reset_index()
        by_config["sem_auc"] = by_config["std_auc"] / by_config["n_permutations"] ** 0.5
        by_config.to_csv(output.by_config, sep="\t", index=False)


rule plot_pos_neg_limit_auc_boxplot:
    input:
        by_config="work_folder/classification/{classifier}/permuted/pos_neg_limit_auc/{network_type}_{esm_model}_auc_by_config.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/pos_neg_limit_auc/plots/{network_type}_{esm_model}_auc_boxplot.png",
    log:
        "logs/classification/{classifier}/permuted/pos_neg_limit_auc/plots/{network_type}_{esm_model}_auc_boxplot.log",
    script:
        "scripts/plot_pos_neg_limit_auc_boxplot.py"


rule compute_permutation_similarity:
    input:
        perm_pos=expand(
            "work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
            permutation=range(config.get("n_permutations", 10)),
        ),
        perm_neg=expand(
            "work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
            permutation=range(config.get("n_permutations", 10)),
        ),
        perm_random=expand(
            "work_folder/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}-random_neg.csv",
            permutation=range(config.get("n_permutations", 10)),
        ),
    output:
        similarity="work_folder/analysis/subset_similarity/permutation_similarity_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}.tsv",
    log:
        "logs/analysis/subset_similarity/permutation_similarity_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}.log",
    script:
        "scripts/compute_permutation_similarity.py"


rule aggregate_permutation_similarity:
    input:
        similarity=expand(
            "work_folder/analysis/subset_similarity/permutation_similarity_{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}.tsv",
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
        ),
    output:
        all_similarity="work_folder/analysis/subset_similarity/all_permutation_similarity_{network_type}.tsv",
    log:
        "logs/analysis/subset_similarity/all_permutation_similarity_{network_type}.log",
    shell:
        """
        echo "dataset\tnetwork_type\tpos_limit\tneg_limit\tlabel\tpermutation_a\tpermutation_b\tjaccard_edges\tspearman_bait\tspearman_prey\tn_a\tn_b" > {output.all_similarity}
        for f in {input.similarity}; do
            tail -n +2 "$f" >> {output.all_similarity}
        done
        """


rule plot_permutation_similarity:
    input:
        similarity="work_folder/analysis/subset_similarity/all_permutation_similarity_{network_type}.tsv",
    output:
        jaccard_plot="work_folder/analysis/subset_similarity/plots/{network_type}_permutation_jaccard.png",
        spearman_plot="work_folder/analysis/subset_similarity/plots/{network_type}_permutation_spearman.png",
    log:
        "logs/analysis/subset_similarity/plots/{network_type}_permutation_similarity.log",
    script:
        "scripts/plot_permutation_similarity.py"


rule evaluate_predictions_random_no_hrni:
    input:
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv",
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_hrni.tsv",
        non_obs="work_folder/subsets/test/{dataset}_{network_type}_selected_no.tsv",
        hrni_saved_model="work_folder/classification/{classifier}/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_model_{esm_model}_parameters.joblib",
        no_saved_model="work_folder/classification/{classifier}/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-random_model_{esm_model}_parameters.joblib",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz"
    output:
        predictions="work_folder/classification/{classifier}/permuted/{permutation}/full_test_predictions/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_model_{esm_model}_predictions.tsv",
    log:
        "logs/classification/{classifier}/permuted/{permutation}/full_test_predictions/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_{esm_model}_predictions.log",
    threads: 1
    resources:
        mem_gb=2,
    script:
        "scripts/evaluate_predictions.py"


rule get_permuted_full_grid_predictions:
    # For one permutation/dataset/network_type, score every trained model in the
    # pos_limit x neg_limit x {train_hrni, train_no} grid against the same three fixed
    # test sets (Interaction, Negative_HRNI, Negative_NO), one prediction column per model.
    # Feeds the future train_hrni-vs-train_no robustness analysis.
    input:
        test_pos="work_folder/subsets/test/{dataset}_{network_type}_pos.csv",
        selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_hrni.tsv",
        no_selected_negative="work_folder/subsets/test/{dataset}_{network_type}_selected_no.tsv",
        models=expand(
            "work_folder/classification/{{classifier}}/permuted/{{permutation}}/model/{{dataset}}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{{esm_model}}_parameters.joblib",
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
            random=["", "-random"],
        ),
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        predictions="work_folder/classification/{classifier}/permuted/{permutation}/full_test_predictions/all_models/{dataset}_{network_type}_model_{esm_model}_predictions.tsv",
    log:
        "logs/classification/{classifier}/permuted/{permutation}/full_test_predictions/all_models/{dataset}_{network_type}_{esm_model}_predictions.log",
    threads: 1
    resources:
        mem_gb=4,
    script:
        "scripts/evaluate_all_model_predictions.py"


rule compute_protein_accuracy:
    input:
        predictions="work_folder/classification/{classifier}/permuted/{permutation}/full_test_predictions/all_models/{dataset}_{network_type}_model_{esm_model}_predictions.tsv",
    output:
        protein_accuracy="work_folder/classification/{classifier}/permuted/{permutation}/full_test_predictions/all_models/protein_accuracy/{dataset}_{network_type}_model_{esm_model}_protein_accuracy.tsv",
    log:
        "logs/classification/{classifier}/permuted/{permutation}/full_test_predictions/all_models/protein_accuracy/{dataset}_{network_type}_{esm_model}_protein_accuracy.log",
    script:
        "scripts/compute_protein_accuracy.py"


rule aggregate_protein_accuracy:
    input:
        protein_accuracy=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/full_test_predictions/all_models/protein_accuracy/{dataset}_{{network_type}}_model_{{esm_model}}_protein_accuracy.tsv",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
        ),
    output:
        all_protein_accuracy="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/all_protein_accuracy_{network_type}_{esm_model}.tsv",
    log:
        "logs/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/all_protein_accuracy_{network_type}_{esm_model}.log",
    run:
        # Pool n_correct/n_total across all 10 permutations (each already summed over that
        # permutation's 6 pos_limit x neg_limit models), so ratio = n_correct / 60 per protein.
        df = pd.concat(
            (pd.read_csv(f, sep="\t", dtype={"protein": "string"}) for f in input.protein_accuracy),
            ignore_index=True,
        )
        count_columns = ["n_correct_interaction", "n_total_interaction", "n_correct_noninteraction", "n_total_noninteraction"]
        summed = df.groupby(["dataset", "network_type", "esm_model", "train_type", "protein"], as_index=False)[count_columns].sum()
        summed["interaction_ratio"] = summed["n_correct_interaction"] / summed["n_total_interaction"]
        summed["noninteraction_ratio"] = summed["n_correct_noninteraction"] / summed["n_total_noninteraction"]
        summed.to_csv(output.all_protein_accuracy, sep="\t", index=False)


rule plot_protein_accuracy_scatter:
    input:
        protein_accuracy="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/all_protein_accuracy_{network_type}_{esm_model}.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/plots/{network_type}_{esm_model}_{dataset}_protein_accuracy_scatter.png",
    log:
        "logs/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/plots/{network_type}_{esm_model}_{dataset}_protein_accuracy_scatter.log",
    script:
        "scripts/plot_protein_accuracy_scatter.py"


rule aggregate_protein_split_robustness:
    # Unlike aggregate_protein_accuracy (which pools n_correct/n_total across all 10 permutations
    # into one ratio per protein), this keeps each permutation's own ratio (already pooled over
    # that permutation's 6 pos_limit x neg_limit configs) and takes the std across the 10
    # permutations per protein/train_type -- i.e. how consistent a protein's classification is
    # specifically across different train/test splits, isolated from pos_limit/neg_limit choice.
    input:
        protein_accuracy=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/full_test_predictions/all_models/protein_accuracy/{dataset}_{{network_type}}_model_{{esm_model}}_protein_accuracy.tsv",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
        ),
    output:
        robustness="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/protein_split_robustness_{network_type}_{esm_model}.tsv",
    log:
        "logs/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/protein_split_robustness_{network_type}_{esm_model}.log",
    run:
        df = pd.concat(
            (pd.read_csv(f, sep="\t", dtype={"protein": "string"}) for f in input.protein_accuracy),
            ignore_index=True,
        )
        df["interaction_ratio"] = df["n_correct_interaction"] / df["n_total_interaction"]
        df["noninteraction_ratio"] = df["n_correct_noninteraction"] / df["n_total_noninteraction"]

        robustness = df.groupby(
            ["dataset", "network_type", "esm_model", "train_type", "protein"], as_index=False
        ).agg(
            interaction_ratio_mean=("interaction_ratio", "mean"),
            interaction_ratio_std=("interaction_ratio", "std"),
            noninteraction_ratio_mean=("noninteraction_ratio", "mean"),
            noninteraction_ratio_std=("noninteraction_ratio", "std"),
            n_permutations=("permutation", "nunique"),
        )
        robustness.to_csv(output.robustness, sep="\t", index=False)


rule plot_protein_split_robustness_boxplot:
    input:
        robustness="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/protein_split_robustness_{network_type}_{esm_model}.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/plots/{network_type}_{esm_model}_protein_split_robustness_boxplot.png",
    log:
        "logs/classification/{classifier}/permuted/full_test_predictions/all_models/protein_accuracy/plots/{network_type}_{esm_model}_protein_split_robustness_boxplot.log",
    script:
        "scripts/plot_protein_split_robustness_boxplot.py"


rule get_prediction_jaccard:
    input:
        predictions = expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/full_test_predictions/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_model_{{esm_model}}_predictions.tsv", permutation=range(config.get("n_permutations", 10))
        ),
    output:
        jaccard = "work_folder/classification/{classifier}/full_test_predictions/jaccard/jaccard_similarity_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_{esm_model}_{pair_set}.tsv",
    log:
        "logs/classification/{classifier}/full_test_predictions/jaccard_similarity_{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_{esm_model}_{pair_set}.log",
    script:
        "scripts/get_prediction_jaccard.py"


rule aggregate_jaccard:
    input:
        jaccard = expand(
            "work_folder/classification/{{classifier}}/full_test_predictions/jaccard/jaccard_similarity_{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}_{{esm_model}}_{{pair_set}}.tsv",
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
        ),
    output:
        all_jaccard = "work_folder/classification/{classifier}/full_test_predictions/jaccard/all_jaccard_{esm_model}_{network_type}_{pair_set}_similarity.tsv",
    log:
        "logs/classification/{classifier}/full_test_predictions/all_jaccard_{esm_model}_{network_type}_{pair_set}_similarity.log",
    shell:
        """
        echo "dataset\tnetwork_type\tpair_set\tpos_limit\tneg_limit\tpermutation\tlabel\tjaccard\tn_pairs" > {output.all_jaccard}
        for f in {input.jaccard}; do
            tail -n +2 "$f" >> {output.all_jaccard}
        done
        """


rule aggregate_negative_accuracy_permuted:
    input:
        hrni=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{{random}}_model_{{esm_model}}_metrics_hrni.txt",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
        ),
        no=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{{random}}_model_{{esm_model}}_metrics_no.txt",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
        ),
    output:
        table="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}{random}_negative_accuracy.tsv",
    log:
        "logs/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}{random}_negative_accuracy.log",
    run:
        def parse_negative_accuracy(path):
            with open(path) as f:
                for line in f:
                    if line.startswith("Accuracy non-interaction:"):
                        return float(line.split(":", 1)[1].strip())
            raise ValueError(f"missing negative accuracy in {path}")

        def parse_key(path, negative_suffix):
            permutation = path.split("/permuted/")[1].split("/")[0]
            stem = path.split("/")[-1].removesuffix(f"_metrics_{negative_suffix}.txt")
            stem = stem.removesuffix(f"_model_{wildcards.esm_model}")
            dataset, rest = stem.split(f"_{wildcards.network_type}_limit_", 1)
            neg_limit, pos_limit = rest.split("_poslim_", 1)
            pos_limit = pos_limit.removesuffix(wildcards.random)
            return dataset, pos_limit, neg_limit, permutation

        no_accuracy = {
            parse_key(path, "no"): parse_negative_accuracy(path) for path in input.no
        }

        with open(output.table, "w") as w:
            w.write("dataset\tpos_limit\tneg_limit\tpermutation\tacc_hrni\tacc_no\n")
            for path in input.hrni:
                key = parse_key(path, "hrni")
                dataset, pos_limit, neg_limit, permutation = key
                acc_hrni = parse_negative_accuracy(path)
                acc_no = no_accuracy[key]
                w.write(f"{dataset}\t{pos_limit}\t{neg_limit}\t{permutation}\t{acc_hrni}\t{acc_no}\n")


rule plot_negative_accuracy_hrni_vs_no:
    input:
        non_random="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}_negative_accuracy.tsv",
        random="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}-random_negative_accuracy.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/negative_accuracy/plots/{network_type}_{esm_model}_negative_accuracy_hrni_vs_no.png",
    log:
        "logs/classification/{classifier}/permuted/negative_accuracy/plots/{network_type}_{esm_model}_negative_accuracy_hrni_vs_no.log",
    script:
        "scripts/plot_negative_accuracy_hrni_vs_no.py"


rule aggregate_hrni_accuracy_permuted:
    input:
        hrni=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{{random}}_model_{{esm_model}}_metrics_hrni.txt",
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
        ),
    output:
        table="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}{random}_hrni_accuracy.tsv",
    log:
        "logs/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}{random}_hrni_accuracy.log",
    run:
        def parse_accuracy(path, label):
            with open(path) as f:
                for line in f:
                    if line.startswith(label):
                        return float(line.split(":", 1)[1].strip())
            raise ValueError(f"missing '{label}' in {path}")

        def parse_key(path):
            permutation = path.split("/permuted/")[1].split("/")[0]
            stem = path.split("/")[-1].removesuffix("_metrics_hrni.txt")
            stem = stem.removesuffix(f"_model_{wildcards.esm_model}")
            dataset, rest = stem.split(f"_{wildcards.network_type}_limit_", 1)
            neg_limit, pos_limit = rest.split("_poslim_", 1)
            pos_limit = pos_limit.removesuffix(wildcards.random)
            return dataset, pos_limit, neg_limit, permutation

        with open(output.table, "w") as w:
            w.write("dataset\tpos_limit\tneg_limit\tpermutation\tacc_pos\tacc_neg\n")
            for path in input.hrni:
                dataset, pos_limit, neg_limit, permutation = parse_key(path)
                acc_pos = parse_accuracy(path, "Accuracy interaction:")
                acc_neg = parse_accuracy(path, "Accuracy non-interaction:")
                w.write(f"{dataset}\t{pos_limit}\t{neg_limit}\t{permutation}\t{acc_pos}\t{acc_neg}\n")


rule plot_hrni_positive_vs_negative_accuracy:
    input:
        non_random="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}_hrni_accuracy.tsv",
        random="work_folder/classification/{classifier}/permuted/negative_accuracy/{network_type}_{esm_model}-random_hrni_accuracy.tsv",
    output:
        plot="work_folder/classification/{classifier}/permuted/negative_accuracy/plots/{network_type}_{esm_model}_hrni_positive_vs_negative_accuracy.png",
    log:
        "logs/classification/{classifier}/permuted/negative_accuracy/plots/{network_type}_{esm_model}_hrni_positive_vs_negative_accuracy.log",
    script:
        "scripts/plot_hrni_positive_vs_negative_accuracy.py"


rule plot_jaccard_negative_hrni_vs_no:
    input:
        jaccard="work_folder/classification/{classifier}/full_test_predictions/jaccard/all_jaccard_{esm_model}_{network_type}_{pair_set}_similarity.tsv",
    output:
        plot="work_folder/classification/{classifier}/full_test_predictions/jaccard/plots/{esm_model}_{network_type}_{pair_set}_jaccard_negative_hrni_vs_no.png",
    log:
        "logs/classification/{classifier}/full_test_predictions/jaccard/plots/{esm_model}_{network_type}_{pair_set}_jaccard_negative_hrni_vs_no.log",
    script:
        "scripts/plot_jaccard_negative_hrni_vs_no.py"


rule plot_jaccard_interaction_vs_negative_hrni:
    input:
        jaccard="work_folder/classification/{classifier}/full_test_predictions/jaccard/all_jaccard_{esm_model}_{network_type}_{pair_set}_similarity.tsv",
    output:
        plot="work_folder/classification/{classifier}/full_test_predictions/jaccard/plots/{esm_model}_{network_type}_{pair_set}_jaccard_interaction_vs_negative_hrni.png",
    log:
        "logs/classification/{classifier}/full_test_predictions/jaccard/plots/{esm_model}_{network_type}_{pair_set}_jaccard_interaction_vs_negative_hrni.log",
    script:
        "scripts/plot_jaccard_interaction_vs_negative_hrni.py"


rule summarize_jaccard_gap_by_dataset:
    input:
        jaccard="work_folder/classification/{classifier}/full_test_predictions/jaccard/all_jaccard_{esm_model}_{network_type}_{pair_set}_similarity.tsv",
    output:
        summary="work_folder/classification/{classifier}/full_test_predictions/jaccard/gap_summary/{esm_model}_{network_type}_{pair_set}_interaction_vs_negative_hrni_gap_by_dataset.tsv",
    log:
        "logs/classification/{classifier}/full_test_predictions/jaccard/gap_summary/{esm_model}_{network_type}_{pair_set}_interaction_vs_negative_hrni_gap_by_dataset.log",
    run:
        df = pd.read_csv(input.jaccard, sep="\t")
        pivot = df.pivot_table(
            index=["dataset", "pos_limit", "neg_limit", "permutation"],
            columns="label",
            values="jaccard",
        ).reset_index()
        pivot["gap"] = pivot["Interaction"] - pivot["Negative_HRNI"]

        summary = pivot.groupby("dataset").agg(
            mean_gap=("gap", "mean"),
            lo=("gap", lambda s: s.quantile(0.025)),
            hi=("gap", lambda s: s.quantile(0.975)),
            n=("gap", "count"),
        ).reset_index()
        summary.to_csv(output.summary, sep="\t", index=False)


rule get_hrni_no_jaccard_pairs:
    # Reuses the already-computed full grid predictions (one column per pos_limit x neg_limit x
    # train_type model, per permutation) to build, per threshold config: 10 HRNI-vs-NO jaccard
    # values (same permutation) and 10 choose 2 = 45 NO-vs-NO jaccard values (across permutations),
    # both restricted to predicted-interaction calls on the Negative_HRNI test set.
    input:
        predictions=expand(
            "work_folder/classification/{{classifier}}/permuted/{permutation}/full_test_predictions/all_models/{{dataset}}_{{network_type}}_model_{{esm_model}}_predictions.tsv",
            permutation=range(config.get("n_permutations", 10)),
        ),
    output:
        jaccard="work_folder/classification/{classifier}/full_test_predictions/all_models/jaccard_pairs/{dataset}_{network_type}_model_{esm_model}_hrni_no_jaccard_pairs.tsv",
    log:
        "logs/classification/{classifier}/full_test_predictions/all_models/jaccard_pairs/{dataset}_{network_type}_{esm_model}_hrni_no_jaccard_pairs.log",
    script:
        "scripts/get_hrni_no_jaccard_pairs.py"


rule plot_hrni_no_jaccard_boxplot:
    # 3 dataset pairs (combined/MS/Y2H) of boxplots on Negative_HRNI prediction overlap:
    # HRNI-vs-NO next to NO-vs-NO, pooling all pos_limit x neg_limit threshold configs per box
    # (no jaccard is ever computed between two different threshold configs).
    input:
        jaccard=expand(
            "work_folder/classification/{{classifier}}/full_test_predictions/all_models/jaccard_pairs/{dataset}_{{network_type}}_model_{{esm_model}}_hrni_no_jaccard_pairs.tsv",
            dataset=config["datasets"],
        ),
    output:
        plot="work_folder/classification/{classifier}/full_test_predictions/all_models/jaccard_pairs/plots/{network_type}_{esm_model}_hrni_no_jaccard_boxplot.png",
    log:
        "logs/classification/{classifier}/full_test_predictions/all_models/jaccard_pairs/plots/{network_type}_{esm_model}_hrni_no_jaccard_boxplot.log",
    script:
        "scripts/plot_hrni_no_jaccard_boxplot.py"


