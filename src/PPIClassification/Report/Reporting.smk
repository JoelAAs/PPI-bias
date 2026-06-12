import re
def get_model_data(wc):
    dataset = wc.dataset
    if wc.dataset == "goldensplit":
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
        "work_folder/subsets/train/{selection}/{data}_pos.csv",
        "work_folder/subsets/train/{selection}/{data}_neg.csv",
        "work_folder/subsets/validation/{selection}/{data}_pos.csv",
        "work_folder/subsets/validation/{selection}/{data}_neg.csv",
        "work_folder/subsets/test/{selection}/{data}_pos.csv",
        "work_folder/subsets/test/{selection}/{data}_neg.csv"
    ]


rule generate_balance_report:
    shadow: "shallow"
    params:
        rmd="src/PPIClassification/Report/sample_balance.rmd",
        main_root = workflow.basedir
    input:
        full_pos="work_folder/subsets/{dataset}_full_{pos_limit}_pos.csv",
        full_neg="work_folder/subsets/{dataset}_full_{neg_limit}_neg.csv",
        train_pos="work_folder/subsets/train/{dataset}_limit_{pos_limit}_{partition_name}_pos.csv",
        train_neg="work_folder/subsets/train/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
        train_balanced_pos="work_folder/subsets/train/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_pos.csv",
        train_balanced_neg="work_folder/subsets/train/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
        validate_pos="work_folder/subsets/validation/{dataset}_limit_{pos_limit}_{partition_name}_pos.csv",
        validate_neg="work_folder/subsets/validation/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
        validate_balanced_pos="work_folder/subsets/validation/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_pos.csv",
        validate_balanced_neg="work_folder/subsets/validation/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
        test_pos="work_folder/subsets/test/{dataset}_limit_{pos_limit}_{partition_name}_pos.csv",
        test_neg="work_folder/subsets/test/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
        test_balanced_pos="work_folder/subsets/test/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_pos.csv",
        test_balanced_neg="work_folder/subsets/test/{balance_method}/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
    output:
        html_report="work_folder/subsets/report/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_{balance_method}.nb.html"
    log:
        "logs/subsets/report/{dataset}_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_{balance_method}.log"
    script:
        "render_balance_report.R"


rule generate_classification_report:
    params:
        rmd="src/PPIClassification/Report/classification_report.rmd",
    input:
        all_models = "work_folder/classification/randomforest/metrics/all_metrics.csv",
        directional_metrics = "work_folder/subsets/degree_balance/all_directional.csv",
        undirectional_metrics = "work_folder/subsets/degree_balance/all_undirectional.csv"
    output:
        html_report="work_folder/subsets/report/{dataset}_{network_type}_{model_configuration}_{partition}_roc_curve.nb.html"
    log:
        "logs/subsets/report/{dataset}_{network_type}_{model_configuration}_{partition}_roc_curve.log"
    script:
        "render_classification_report.R"
