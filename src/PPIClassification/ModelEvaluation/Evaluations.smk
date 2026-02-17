
def get_model_validation_data(wc):
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
        f"work_folder{pn}/subsets/validation/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_neg.csv"
    ]

rule get_model_metrics:
    params:
        script_location = "src/PPIClassification/ModelEvaluation/evaluate_model.py"
    input:
        validation_data = lambda wc: get_model_validation_data(wc),
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{model_configuration}}_model_parameters.joblib",
        dummy_baseline = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{model_configuration}}_dummy_baseline_parameters.joblib",
        protein_embeddings = f"work_folder{pn}/embeddings/{{dataset}}_protein_embeddings.csv"
    output:
        metrics=f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_{{model_configuration}}_metrics.txt"
    shell:
        """
        python3 {params.script_location} \
            --validation_ppi_data_pos {input.validation_data[0]} \
            --validation_ppi_data_neg {input.validation_data[1]} \
            --protein_embeddings_file {input.protein_embeddings} \
            --model_file {input.saved_model} \
            --dummy_baseline_file {input.dummy_baseline} \
            --saved_model {input.saved_model} \
            --output_metrics {output.metrics} \
            --threads {threads}
        """