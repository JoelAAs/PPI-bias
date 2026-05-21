rule get_all_balance_metrics:
    input:
        metrics=expand(
            "work_folder{pn}/subsets/train/equal_edge/balance/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_degree.csv",
            pn=pn,
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
            random=["", "-random"],
        ),
    output:
        all_models=f"work_folder{pn}/subsets/train/equal_edge/balance/all_metrics.csv",
    shell:
        """
        echo "dataset\tpos_limit\tneg_limit\trandom\tbait_degree_delta\tprey_degree_delta\tnum_edges" > {output.all_models}
        cat {input.metrics} >> {output.all_models}
        """


rule get_train_degree_balance:
    input:
        train_pos=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
    output:
        balance=f"work_folder{pn}/subsets/train/equal_edge/balance/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_degree.csv",
    run:
        pos_df = pd.read_csv(input.train_pos, sep="\t")
        neg_df = pd.read_csv(input.train_neg, sep="\t")

        pos_G = nx.from_pandas_edgelist(
            pos_df, "bait", "prey", create_using=nx.DiGraph
        )
        neg_G = nx.from_pandas_edgelist(
            neg_df, "bait", "prey", create_using=nx.DiGraph
        )

        pos_out = dict(pos_G.out_degree())
        neg_out = dict(neg_G.out_degree())
        pos_in = dict(pos_G.in_degree())
        neg_in = dict(neg_G.in_degree())

        all_bait = set(pos_out) | set(neg_out)
        all_prey = set(pos_in) | set(neg_in)
        degree_bait_delta = sum(abs(pos_out[n] - neg_out[n]) for n in all_bait)
        degree_prey_delta = sum(abs(pos_in[n] - neg_in[n]) for n in all_prey)
        n_edges = pos_G.number_of_edges() + neg_G.number_of_edges()
        with open(output.balance, "w") as w:
            w.write(
                f"{wildcards.dataset}\t{wildcards.pos_limit}\t{wildcards.neg_limit}\t{wildcards.random!= ""}\t{degree_bait_delta}\t{degree_prey_delta}\t{n_edges}\n"
            )


rule get_model_metrics:
    input:
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_neg.csv",
        saved_model=f"work_folder{pn}/classification/{{classifier}}/model/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        metrics=f"work_folder{pn}/classification/{{classifier}}/metrics/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_metrics.txt",
        roc_png=f"work_folder{pn}/classification/{{classifier}}/metrics/plot/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_roc_curve.png",
    log:
        f"logs{pn}/classification/{{classifier}}/metrics/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_metrics.log",
    threads: 10
    resources:
        mem_gb=20,
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
            --plot_roc_png {output.roc_png} > {log} 2>&1
        """


rule all_metrics:
    input:
        metrics=expand(
            "work_folder{pn}/classification/{{classifier}}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{{esm_model}}_metrics.txt",
            pn=pn,
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
            random=["", "-random"],
        ),
    output:
        all_models=f"work_folder{pn}/classification/{{classifier}}/metrics/all_metrics_{{network_type}}_{{esm_model}}.csv",
    log:
        f"logs{pn}/classification/{{classifier}}/metrics/all_metrics_{{network_type}}_{{esm_model}}.log",
    run:
        with open(output[0], "a") as w:
            w.write("model\troc_auc\tsamples\n")
            for metric_file in input.metrics:
                with open(metric_file, "r") as f:
                    line_out = [line.strip().split(": ")[1] for line in f]
                    line_out = "\t".join(line_out)
                    model_name = metric_file.split("/")[-1].replace(
                        "_metrics.txt", ""
                    )
                    line_out = model_name + "\t" + line_out + "\n"
                    w.write(line_out)


rule get_model_metrics_permuted:
    input:
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_neg.csv",
        saved_model=f"work_folder{pn}/classification/{{classifier}}/permuted/{{permutation}}/model/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        metrics=f"work_folder{pn}/classification/{{classifier}}/permuted/{{permutation}}/metrics/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_metrics.txt",
        roc_png=f"work_folder{pn}/classification/{{classifier}}/permuted/{{permutation}}/metrics/plot/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_roc_curve.png",
    log:
        f"logs{pn}/classification/{{classifier}}/permuted/{{permutation}}/metrics/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_metrics.log",
    threads: 10
    resources:
        mem_gb=80,
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
            --plot_roc_png {output.roc_png} > {log} 2>&1
        """


rule all_metrics_permuted:
    input:
        metrics=expand(
            "work_folder{pn}/classification/{{classifier}}/permuted/{permutation}/metrics/{dataset}_{{network_type}}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{{esm_model}}_metrics.txt",
            pn=pn,
            permutation=range(config.get("n_permutations", 10)),
            dataset=config["datasets"],
            pos_limit=config["positive_limits"],
            neg_limit=config["negative_limits"],
            random=[""],
        ),
    output:
        all_models=f"work_folder{pn}/classification/{{classifier}}/permuted/all_metrics_{{network_type}}_{{esm_model}}.csv",
    log:
        f"logs{pn}/classification/{{classifier}}/permuted/all_metrics_{{network_type}}_{{esm_model}}.log",
    run:
        with open(output[0], "w") as w:
            w.write("permutation\tmodel\troc_auc\tsamples\n")
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


rule compute_train_similarity:
    input:
        similarity_tsv=f"work_folder{pn}/protein_sequences/similarity/sequencesimilarity.tsv",
        perm_pos=expand(
            "work_folder{pn}/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
            pn=pn,
            permutation=range(config.get("n_permutations", 10)),
            ),
        perm_neg=expand(
            "work_folder{pn}/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
            pn=pn,
            permutation=range(config.get("n_permutations", 10))
        ),
        perm_random=expand(
            "work_folder{pn}/subsets/train/permuted/{permutation}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}-random_neg.csv",
            pn=pn,
            permutation=range(config.get("n_permutations", 10))
        )
    output:
        per_protein_similarity=f"work_folder{pn}/analysis/train_similarity/protein_avg_similarity_{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.tsv",
        per_edge_similarity=f"work_folder{pn}/analysis/train_similarity/pair_similarity_{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.tsv"
    log:
        f"logs{pn}/analysis/train_similarity/compute_similarity_{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.log",
    script:
        "compute_train_similarity.py"
        

rule plot_train_similarity:
    input:
        per_protein_similarity=f"work_folder{pn}/analysis/train_similarity/protein_avg_similarity_{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.tsv",
        per_edge_similarity=f"work_folder{pn}/analysis/train_similarity/pair_similarity_{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}.tsv"
    output:
        avg_protein_similarity=f"work_folder{pn}/analysis/train_similarity/plots/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_avg_protein_similarity.png",
        edge_similarity_density=f"work_folder{pn}/analysis/train_similarity/plots/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pair_similarity_density.png",
    log:
        f"logs{pn}/analysis/train_similarity/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_plot.log",
    script: "plot_train_similarity.R"