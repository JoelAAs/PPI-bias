rule random_forest:
    input:
        train_pos=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        params=f"work_folder{pn}/classification/randomforest/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.txt",
        saved_model=f"work_folder{pn}/classification/randomforest/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
    log:
        f"logs{pn}/classification/randomforest/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_model.log"
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
            > {log} 2>&1
        """


rule random_forest_permuted:
    input:
        train_pos=f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        params=f"work_folder{pn}/classification/randomforest/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.txt",
        saved_model=f"work_folder{pn}/classification/randomforest/permuted/{{permutation}}/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
    log:
        f"logs{pn}/classification/randomforest/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_model.log"
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
            > {log} 2>&1
        """


rule xgboost_permuted:
    input:
        train_pos=f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        params=f"work_folder{pn}/classification/xgboost/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.txt",
        saved_model=f"work_folder{pn}/classification/xgboost/permuted/{{permutation}}/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
    log:
        f"logs{pn}/classification/xgboost/permuted/{{permutation}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_model.log"
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
            > {log} 2>&1
        """


rule xgboost:
    input:
        train_pos=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    output:
        params=f"work_folder{pn}/classification/xgboost/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.txt",
        saved_model=f"work_folder{pn}/classification/xgboost/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_{{esm_model}}_parameters.joblib",
    log:
        f"logs{pn}/classification/xgboost/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_{{esm_model}}_model.log"
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
            > {log} 2>&1
        """
