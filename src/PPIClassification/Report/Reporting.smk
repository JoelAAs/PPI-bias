import re
def get_model_data(wc):
    dataset = wc.dataset
    is_mixed = re.match(r"mixed-(.+)_", dataset)
    if is_mixed
        negative_dataset = is_mixed.group(1)
        neg_selection = "undirectionalbalanced"



    return [
        f"work_folder{pn}/subsets/train/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/train/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_neg.csv"
    ]



    elif wc.dataset == "goldensplit":
        data = f"data_{wc.network_type}"
        selection = wc.dataset
    else:
        pos_limit = config["models"][wc.model_configuration]["pos"]
        neg_limit = config["models"][wc.model_configuration]["neg"]
        data =  f"{wc.dataset}_{wc.network_type}_limit_{neg_limit}_poslim_{pos_limit}_{wc.partition}"
        if wc.network_type == "directional":
            selection = "maxflow"
        elif wc.network_type == "undirectional":
            selection="undirectionalbalanced"
        else:
            raise ValueError(f"unkown networktype {wc.network_type}")
    return [
        f"work_folder{pn}/subsets/train/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/train/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/validation/{selection}/{data}_neg.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_pos.csv",
        f"work_folder{pn}/subsets/test/{selection}/{data}_neg.csv"
    ]


rule generate_balance_report:
    shadow: "shallow"
    params:
        rmd="src/PPIClassification/Report/sample_balance.rmd",
        main_root = workflow.basedir
    input:
        full_pos=f"work_folder{pn}/subsets/{{dataset}}_full_{{pos_limit}}_pos.csv",
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_full_{{neg_limit}}_neg.csv",
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        train_balanced_pos=f"work_folder{pn}/subsets/train/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        train_balanced_neg=f"work_folder{pn}/subsets/train/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        validate_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        validate_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        validate_balanced_pos=f"work_folder{pn}/subsets/validation/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        validate_balanced_neg=f"work_folder{pn}/subsets/validation/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        test_balanced_pos=f"work_folder{pn}/subsets/test/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        test_balanced_neg=f"work_folder{pn}/subsets/test/{{balance_method}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
    output:
        html_report=f"work_folder{pn}/subsets/report/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_{{balance_method}}.nb.html"
    script:
        "render_balance_report.R"


rule generate_classification_report:
    shadow: "shallow"
    params:
        rmd="src/PPIClassification/Report/classification_report.rmd",
        main_root = workflow.basedir
    input:
        data = lambda wc: get_model_data(wc),
        pr_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_pr_curve.png",
        pr_neg_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_pr_neg_curve.png",
        roc_png=f"work_folder{pn}/classification/randomforest/metrics/plot/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_roc_curve.png"

    output:
        html_report=f"work_folder{pn}/subsets/report/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_roc_curve.nb.html"
    script:
        "render_classification_report.R"
