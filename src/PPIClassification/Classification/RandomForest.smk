rule random_forest:
    input:
        train_pos=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/equal_edge/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        validation_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_pos.csv",
        validation_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_neg.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz",
    output:
        params=f"work_folder{pn}/classification/randomforest/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.txt",
        saved_model=f"work_folder{pn}/classification/randomforest/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.joblib",
    threads: 23
    params:
        script_location="src/PPIClassification/Classification/ppi_classify_rf.py",
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
