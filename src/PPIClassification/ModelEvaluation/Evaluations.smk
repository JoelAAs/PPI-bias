def input_metrics(wc):
    expected_input = expand(
        f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_{{model_configuration}}_{{partition}}_metrics.txt",
         dataset=config["datasets"], model_configuration=(c for c in config["models"].keys() if c != "goldensplit"))
    expected_input.append(f"work_folder{pn}/classification/randomforest/metrics/goldensplit_asis_metrics.txt")
    return expected_input

def get_model_validation_data(wc):
    if wc.dataset == "goldensplit":
        data = "data"
        selection = wc.dataset
    else:
        pos_limit = config["models"][wc.model_configuration]["pos"]
        neg_limit = config["models"][wc.model_configuration]["neg"]
        data =  f"{wc.dataset}_limit_{neg_limit}_poslim_{pos_limit}_{wc.partition}"
        if wc.network_type == "directional":
            selection = "maxflow"
        elif wc.network_type == "undirectional":
            selection="undirectionalbalanced"
    return [
        f"work_folder{pn}/subsets/test/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_neg.csv"
    ]

rule get_model_metrics:
    params:
        script_location = "src/PPIClassification/ModelEvaluation/evaluate_model.py"
    input:
        validation_data = lambda wc: get_model_validation_data(wc),
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{model_configuration}}_{{partition}}_model_parameters.joblib",
        protein_embeddings = f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        metrics=f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_{{model_configuration}}_{{partition}}_metrics.txt",
        pr_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{model_configuration}}_{{partition}}_pr_curve.png",
        pr_neg_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{model_configuration}}_{{partition}}_pr_neg_curve.png",
        roc_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{model_configuration}}_{{partition}}_roc_curve.png"
    shell:
        """
        python3 {params.script_location} \
            --pos_data_file {input.validation_data[0]} \
            --neg_data_file {input.validation_data[1]} \
            --protein_embeddings_file {input.protein_embeddings} \
            --model_file {input.saved_model} \
            --output_file {output.metrics} \
            --plot_pr_png {output.pr_png} \
            --plot_neg_pr_png {output.pr_neg_png} \
            --plot_roc_png {output.roc_png}
        """

rule all_metrics:
    input:
        metrics = lambda wc: input_metrics(wc)    
    output:
        expected_input = expand(
            f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_{{model_configuration}}_{{partition}}_metrics.txt",
         dataset=config["datasets"], model_configuration=config["models"], partition=config["partitions"])
        run:
        with open(output[0], "a") as w:
            w.write("model\tpr_auc\tpr_auc_base\tpr_auc_neg\tpr_auc_neg_base\troc_auc\troc_auc_base\n")
            for metric_file in input.metrics:
                with open(metric_file, "r") as f:
                    line_out = [line.strip().split(": ")[1] for line in f]
                    line_out = "\t".join(line_out) 
                    model_name = metric_file.split("/")[-1].replace("_metrics.txt", "")
                    line_out = model_name + "\t" + line_out + "\n"
                    w.write(line_out)


