

rule get_model_metrics:
    params:
        script_location = "src/PPIClassification/ModelEvaluation/evaluate_model.py"
    input:
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_directional_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_directional_neg.csv",
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_model_parameters.joblib",
        protein_embeddings = f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    output:
        metrics=f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_metrics.txt",
        pr_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_pr_curve.png",
        pr_neg_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_pr_neg_curve.png",
        ce_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_ce.png",
    threads: 10
    log:
        f"logs{pn}/classification/randomforest/metrics/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_metrics.log"
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
            --plot_ce_png {output.ce_png} > {log} 2>&1
        """

rule all_metrics:
    input:
        metrics = expand(
            f"work_folder{pn}/classification/randomforest/metrics/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}{{random}}_metrics.txt",
            dataset=config["datasets"], pos_limit=config["positive_limits"], neg_limit=config["negative_limits"], random=["", "-random"])
    output:
        all_models = f"work_folder{pn}/classification/randomforest/metrics/all_metrics.csv"
    log:
        f"logs{pn}/classification/randomforest/metrics/all_metrics.log"
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

