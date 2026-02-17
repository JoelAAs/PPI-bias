
def get_expected_input(wc):
    if wc.dataset in config["models"]:
        data = "data"
        selection = wc.dataset
    else:
        pos_limit = config["models"][wc.model_configuration]["pos"]
        neg_limit = config["models"][wc.model_configuration]["neg"]
        selection = config["models"][wc.model_configuration]["balancing"]
        partition_name = config["models"][wc.model_configuration]["partition"]
        data =  f"{wc.dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}"

    return [
        f"work_folder{pn}/subsets/train/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/train/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_neg.csv"
    ]


rule random_forest:
    params:
        script_location = "src/PPIClassification/Classification/ppi_classify_rf.py"
    input:
        data = lambda wc: get_expected_input(wc),
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        params=f"work_folder{pn}/classification/randomforest/{{dataset}}_{{model_configuration}}_model_parameters.txt",
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{model_configuration}}_model_parameters.joblib",
        baseline_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{model_configuration}}_dummy_model.joblib"
    threads: 48
    shell:
        """
        python3 {params.script_location} \
            --train_ppi_data_pos {input.data[0]} \
            --train_ppi_data_neg {input.data[1]} \
            --validation_ppi_data_pos {input.data[2]} \
            --validation_ppi_data_neg {input.data[3]} \
            --test_ppi_data_pos {input.data[4]} \
            --test_ppi_data_neg {input.data[5]} \
            --protein_embeddings {input.protein_embeddings} \
            --params_out {output.params} \
            --threads {threads} \
            --randomstate 1234 \
            --saved_model {output.saved_model} \
            --saved_dummy_classifer {output.baseline_model}
        """