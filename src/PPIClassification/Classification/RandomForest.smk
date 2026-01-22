rule random_forest_max_flow:
    params:
        script_location = "src/PPIClassification/ppi_classify_rf"
    input:
        train_ppi_data_pos=f"work_folder{pn}/subsets/train/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        train_ppi_data_neg=f"work_folder{pn}/subsets/train/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        validation_ppi_data_pos=f"work_folder{pn}/subsets/validation/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        validation_ppi_data_neg=f"work_folder{pn}/subsets/validation/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        test_ppi_data_pos=f"work_folder{pn}/subsets/test/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_ppi_data_neg=f"work_folder{pn}/subsets/test/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        params= f"work_folder{pn}/classification/randomforest/parametes_{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.txt"
    threads: 48
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.train_ppi_data_pos} \
            --train_ppi_data_neg {input.train_ppi_data_neg} \
            --validation_ppi_data_pos {input.validation_ppi_data_pos} \
            --validation_ppi_data_neg {input.validation_ppi_data_neg} \
            --test_ppi_data_pos {input.test_ppi_data_pos} \
            --test_ppi_data_neg {input.test_ppi_data_neg} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234
        """