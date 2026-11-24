rule naive_xgboost:
    input:
        train="work_folder/naive_split/seed_{seed}_train/{network_type}/balanced/{pos_set}_{neg_set}{no}.tsv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/naive_split/xgboost/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/naive_split/xgboost/model/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}.joblib",
    log:
        "logs/classification/naive_split/xgboost/{network_type}_{pos_set}_{neg_set}{no}_seed_{seed}_{esm_model}_model.log"
    threads: 10
    resources:
        mem_gb=40
    params:
        random_state=lambda wildcards: int(wildcards.seed),
        predef=False
    script:
        "scripts/ppi_classify_xgboost.py"

rule predef_xgboost:
    input:
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        protein_embeddings = "work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz"
    output:
        params="work_folder/classification/predef/xgboost/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/predef/xgboost/model/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
    log:
        "logs/classification/predef/xgboost/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_model.log",
    threads: 10
    resources:
        mem_gb=40,
    params:
        random_state=1234,
        predef=True
    script:
        "scripts/ppi_classify_xgboost.py"

