rule random_forest:
    input:
        train_pos="work_folder/subsets/train/equal_edge/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/train/equal_edge/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
        validation_pos="work_folder/subsets/validation/{dataset}_{network_type}_pos.csv",
        validation_neg="work_folder/subsets/validation/{dataset}_{network_type}_neg.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/randomforest/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/randomforest/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
    log:
        "logs/classification/randomforest/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_model.log"
    threads: 15
    resources:
        mem_gb=40
    params:
        script_location="src/PPIClassification/Classification/ppi_classify_rf.py",
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_pos} \
            --train_ppi_data_neg {input.train_neg} \
            --validation_ppi_data_pos {input.validation_pos} \
            --validation_ppi_data_neg {input.validation_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} \
            --network_type {wildcards.network_type} \
            > {log} 2>&1
        """


rule random_forest_permuted:
    input:
        train_pos="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
        validation_pos="work_folder/subsets/validation/{dataset}_{network_type}_pos.csv",
        validation_neg="work_folder/subsets/validation/{dataset}_{network_type}_neg.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/randomforest/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/randomforest/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
    log:
        "logs/classification/randomforest/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_model.log"
    threads: 15
    resources:
        mem_gb=40
    params:
        script_location="src/PPIClassification/Classification/ppi_classify_rf.py",
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_pos} \
            --train_ppi_data_neg {input.train_neg} \
            --validation_ppi_data_pos {input.validation_pos} \
            --validation_ppi_data_neg {input.validation_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} \
            --network_type {wildcards.network_type} \
            > {log} 2>&1
        """


rule xgboost_permuted:
    input:
        train_pos="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/train/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
        validation_pos="work_folder/subsets/validation/{dataset}_{network_type}_pos.csv",
        validation_neg="work_folder/subsets/validation/{dataset}_{network_type}_neg.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/xgboost/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/xgboost/permuted/{permutation}/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
    log:
        "logs/classification/xgboost/permuted/{permutation}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_model.log"
    threads: 15
    resources:
        mem_gb=40
    params:
        script_location="src/PPIClassification/Classification/ppi_classify_xgboost.py",
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_pos} \
            --train_ppi_data_neg {input.train_neg} \
            --validation_ppi_data_pos {input.validation_pos} \
            --validation_ppi_data_neg {input.validation_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} \
            --network_type {wildcards.network_type} \
            > {log} 2>&1
        """


rule xgboost:
    input:
        train_pos="work_folder/subsets/train/equal_edge/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        train_neg="work_folder/subsets/train/equal_edge/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_neg.csv",
        validation_pos="work_folder/subsets/validation/{dataset}_{network_type}_pos.csv",
        validation_neg="work_folder/subsets/validation/{dataset}_{network_type}_neg.csv",
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    output:
        params="work_folder/classification/xgboost/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/xgboost/model/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_model_{esm_model}_parameters.joblib",
    log:
        "logs/classification/xgboost/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}{random}_{esm_model}_model.log"
    threads: 15
    resources:
        mem_gb=40
    params:
        script_location="src/PPIClassification/Classification/ppi_classify_xgboost.py",
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_pos} \
            --train_ppi_data_neg {input.train_neg} \
            --validation_ppi_data_pos {input.validation_pos} \
            --validation_ppi_data_neg {input.validation_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} \
            --network_type {wildcards.network_type} \
            > {log} 2>&1
        """
