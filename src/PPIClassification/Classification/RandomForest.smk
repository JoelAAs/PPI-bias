
rule random_forest:
    params:
        script_location = "src/PPIClassification/Classification/ppi_classify_rf.py"
    input:
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        params =      f"work_folder{pn}/classification/randomforest/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.txt",
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.joblib"
    threads: 23
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_pos} \
            --train_ppi_data_neg {input.train_neg} \
            --validation_ppi_data_pos {input.validation_pos} \
            --validation_ppi_data_neg {input.validation_neg} \
            --test_ppi_data_pos {input.test_pos} \
            --test_ppi_data_neg {input.test_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} 
        """



rule get_degree_balance:
    params:
        script_location = "src/PPIClassification/ModelEvaluation/degree_balance_metrics.py"
    input:    
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv"
    output:
        degree_balance = f"work_folder{pn}/subsets/degree_balance/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}.csv"
    shell:
        """
        python3 {params.script_location} \
            --pos_train {input.train_pos} \
            --neg_train {input.train_neg} \
            --pos_val {input.validation_pos} \
            --neg_val {input.validation_neg} \
            --pos_test {input.test_pos} \
            --neg_test  {input.test_neg} \
            --output_file {output.degree_balance} \
            --network_type {wildcards.network_type}
        """
