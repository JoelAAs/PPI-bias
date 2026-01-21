rule random_forest:
    input:
        train_ppi_data_pos=f"work_folder{pn}/subsets/train/balanced/{{dataset}}_pos.csv",
        train_ppi_data_neg=f"work_folder{pn}/subsets/train/balanced/{{dataset}}_neg.csv",
        test_ppi_data_pos=f"work_folder{pn}/subsets/test/balanced/{{dataset}}_pos.csv",
        test_ppi_data_neg=f"work_folder{pn}/subsets/test/balanced/{{dataset}}_neg.csv",
        validation_ppi_data_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_pos.csv",
        validation_ppi_data_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_neg.csv",
        protein_embedings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        params= f"work_folder{pn}/classification/randomforest/parametes_{{dataset}}.txt"
