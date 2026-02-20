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