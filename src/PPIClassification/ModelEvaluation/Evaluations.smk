def input_metrics():
    metrics = expand(
        f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_metrics.txt",
        dataset=config["datasets"], neg_limit=1,pos_limit=[0.02, 0.15, 0.29], random=["", "-random"])

    metrics = [m for m in metrics if m not in failed]
    return matrics


rule get_model_metrics:
    params:
        script_location = "src/PPIClassification/ModelEvaluation/evaluate_model.py"
    input:
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_neg.csv",
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.joblib",
        protein_embeddings = f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        metrics=f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_metrics.txt",
        pr_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_pr_curve.png",
        pr_neg_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_pr_neg_curve.png",
        roc_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_roc_curve.png",
        ce_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_ce.png",
    shell:
        """
        python3 {params.script_location} \
            --pos_data_file {input.test_pos} \
            --neg_data_file {input.test_neg} \
            --protein_embeddings_file {input.protein_embeddings} \
            --model_file {input.saved_model} \
            --output_file {output.metrics} \
            --plot_pr_png {output.pr_png} \
            --plot_neg_pr_png {output.pr_neg_png} \
            --plot_roc_png {output.roc_png} \
            --plot_ce_png {output.ce_png}
        """

rule all_metrics:
    input:
        input_metrics
    output:
        all_models = f"work_folder{pn}/classification/randomforest/metrics/all_metrics.csv"
    run:
        with open(output[0], "a") as w:
            w.write("model\tpr_auc\tpr_auc_base\tpr_auc_neg\tpr_auc_neg_base\troc_auc\troc_auc_base\tce_obs\tce_baseline\tsamples\n")
            for metric_file in input.metrics:
                with open(metric_file, "r") as f:
                    line_out = [line.strip().split(": ")[1] for line in f]
                    line_out = "\t".join(line_out) 
                    model_name = metric_file.split("/")[-1].replace("_metrics.txt", "")
                    line_out = model_name + "\t" + line_out + "\n"
                    w.write(line_out)

