

rule predef_xgboost:
    input:
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/predef/xgboost/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/predef/xgboost/model/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
    log:
        "logs/classification/predef/xgboost/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_model.log",
    threads: 10
    resources:
        mem_gb=40,
    params:
        random_state=1234,
    script:
        "scripts/ppi_classify_xgboost_naive.py"


rule predef_model_metrics:
    input:
        saved_model="work_folder/classification/predef/xgboost/model/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        test="results_matched/{pos_set}_{neg_set}-{test_state}/test_balanced.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        metrics="work_folder/classification/predef/xgboost/metrics/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
        pred="work_folder/classification/predef/xgboost/predictions/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_pred.tsv.gz",
    log:
        "logs/classification/predef/xgboost/metrics/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.log",
    threads: 1
    resources:
        mem_gb=8,
    script:
        "../ModelEvaluation/scripts/evaluate_predefined_split.py"


rule predef_all_metrics:
    input:
        metrics=lambda wc: expand(
            "work_folder/classification/predef/xgboost/metrics/undirectional_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
            esm_model = 
            pos_set=["hippie", "fdr", "pod"],
            neg_set=["fdr", "pod"],
            state=["cand", "nonobs"],
            test_state=["cand", "nonobs"]
        ),
    output:
        all_metrics="work_folder/classification/predef/xgboost/all_metrics_undirectional.tsv",
    shell:
        "awk 'FNR > 1 || NR == 1' {input.metrics} > {output.all_metrics}"
