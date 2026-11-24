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


## Naive split
rule naive_model_metrics:
    params:
        predef=False
    input:
        saved_model="work_folder/classification/naive_split/{classifier}/model/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}.joblib",
        train="work_folder/naive_split/seed_{seed}_train/{network_type}/balanced/{pos_set}_{neg_set}{no}.tsv",
        test="work_folder/naive_split/seed_{seed}_test/{network_type}/balanced/{pos_set}_{neg_set}{test_no}.tsv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        metrics="work_folder/classification/naive_split/{classifier}/metrics/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_test{test_no}_metrics.tsv",
        pred="work_folder/classification/naive_split/{classifier}/predictions/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_test{test_no}_pred.tsv.gz",
    log:
        "logs/classification/naive_split/{classifier}/metrics/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_test{test_no}_metrics.log",
    threads: 1
    resources:
        mem_gb=3,
    script:
        "scripts/evaluate_classification.py"


rule naive_all_metrics:
    input:
        metrics=lambda wc: expand(
            "work_folder/classification/naive_split/{classifier}/metrics/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_test{test_no}_metrics.tsv",
            **dict(wc.items()),  # classifier, network_type, esm_model
            pos_set=["hippie"] + list(config["naive_split"]["positive"]),
            neg_set=list(config["naive_split"]["negative"]),
            no=["", "-no"],
            test_no=["", "-no"],
            seed=config["naive_split"]["seeds"],
        ),
    output:
        all_metrics="work_folder/classification/naive_split/{classifier}/all_metrics_{network_type}_{esm_model}.tsv",
    shell:
        "awk 'FNR > 1 || NR == 1' {input.metrics} > {output.all_metrics}"

rule plot_naive_performance:
    input:
        all_metrics=expand("work_folder/classification/naive_split/{{classifier}}/all_metrics_{{network_type}}_{esm_model}.tsv",
            esm_model=config["naive_split"]["esm_model"]
        )
    output:
        auc="work_folder/classification/naive_split/plots/{network_type}_{classifier}_auc_hrni_test.svg",
        auc_embedding="work_folder/classification/naive_split/plots/{network_type}_{classifier}_auc_hrni_test_by_embedding.svg",
        delta_test="work_folder/classification/naive_split/plots/{network_type}_{classifier}_delta_hrni_vs_no_test.svg",
        delta_adj="work_folder/classification/naive_split/plots/{network_type}_{classifier}_delta_auc_vs_adjusted.svg",
    script:
        "scripts/plot_naive_performance.R"

## Predef
rule predef_model_metrics:
    params:
        predef=True
    input:
        saved_model="work_folder/classification/predef/{classifier}/model/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        test="results_matched/{pos_set}_{neg_set}-{test_state}/test_balanced.csv",
        protein_embeddings = "work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz"
    output:
        metrics="work_folder/classification/predef/{classifier}/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
        pred="work_folder/classification/predef/{classifier}/predictions/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_pred.tsv.gz",
    log:
        "logs/classification/predef/{classifier}/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.log",
    threads: 1
    resources:
        mem_gb=8,
    script:
        "scripts/evaluate_classification.py"


rule predef_all_metrics:
    input:
        metrics=lambda wc: expand(
            "work_folder/classification/predef/xgboost/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
            network_type="undirectional",
            esm_model = config["naive_split"]["esm_model"],
            pos_set=["hippie", "fdr", "pod"],
            neg_set=["fdr", "pod"],
            state=["cand", "nonobs"],
            test_state=["cand", "nonobs"]
        ),
    output:
        all_metrics="work_folder/classification/predef/xgboost/all_metrics.tsv",
    shell:
        "awk 'FNR > 1 || NR == 1' {input.metrics} > {output.all_metrics}"


# rule plot_pair_embedding_cosine_similarity:
#     params:
#         score_col=lambda wildcards: {"fdr": "log_odds", "pod": "lower_bound_pod"}[wildcards.source],
#         bait_col=f"{config['id_pattern']}_bait",
#         prey_col=f"{config['id_pattern']}_prey",
#         neg_limit=lambda wildcards:config["classification_sources"][wildcards.source]["negative"]["limits"][0],
#     input:
#         dataset=lambda wildcards: config["classification_sources"][wildcards.source]["path"].format(
#             network_type=wildcards.network_type, dataset=wildcards.dataset),
#         protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
#     output:
#         emb_sim_vs_score="work_folder/classification/embedding_sim_auc/plots/{dataset}_{source}_{network_type}_{esm_model}_cosine_sim_vs_score.png",
#     log:
#         "logs/classification/{source}/embedding_sim_auc/plots/{dataset}_{network_type}_{esm_model}_cosine_sim_vs_score.log",
#     resources:
#         mem_gb=16,
#     script:
#         "scripts/plot_cosine_sim_vs_score.py"


rule compare_predef_naive:
    input:
        predef_metrics="work_folder/classification/predef/xgboost/all_metrics.tsv",
        naive_metrics= expand("work_folder/classification/naive_split/xgboost/all_metrics_undirectional_{esm_model}.tsv",
            esm_model=config["naive_split"]["esm_model"]
        )
    output:
        pr_na="work_folder/classification/performance/predef_and_naive.tsv"
    shell:
        "awk 'FNR > 1 || NR == 1' {input} > {output}"
