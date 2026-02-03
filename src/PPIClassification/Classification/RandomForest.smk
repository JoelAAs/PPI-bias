rule random_forest:
    params:
        script_location = "src/PPIClassification/Classification/ppi_classify_rf.py"
    input:
        train_ppi_data_pos=f"work_folder{pn}/subsets/train/{{selection}}/{{dataset}}_{{parameters}}_pos.edgelist",
        train_ppi_data_neg=f"work_folder{pn}/subsets/train/{{selection}}/{{dataset}}_{{parameters}}_neg.edgelist",
        validation_ppi_data_pos=f"work_folder{pn}/subsets/validation/{{selection}}/{{dataset}}_{{parameters}}_pos.edgelist",
        validation_ppi_data_neg=f"work_folder{pn}/subsets/validation/{{selection}}/{{dataset}}_{{parameters}}_neg.edgelist",
        test_ppi_data_pos=f"work_folder{pn}/subsets/test/{{selection}}/{{dataset}}_{{parameters}}_pos.edgelist",
        test_ppi_data_neg=f"work_folder{pn}/subsets/test/{{selection}}/{{dataset}}_{{parameters}}_neg.edgelist",
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        params= f"work_folder{pn}/classification/randomforest/{{selection}}/{{dataset}}_{{parameters}}_model_parameters.txt"
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