rule generate_balance_report:
    params:
        rmd_location = "src/PPIClassification/sample_balance.rmd"
    input:
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_full_{{pos_limit}}_pos.csv",
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_full_{{neg_limit}}_neg.csv",
        train_pos = f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        train_neg =f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        train_maxflow_pos= f"work_folder{pn}/subsets/train/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        train_maxflow_neg= f"work_folder{pn}/subsets/train/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        validate_pos = f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        validate_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        validate_maxflow_pos= f"work_folder{pn}/subsets/validation/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        validate_maxflow_neg= f"work_folder{pn}/subsets/validation/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        test_pos = f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv",
        test_maxflow_pos=f"work_folder{pn}/subsets/test/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        test_maxflow_neg=f"work_folder{pn}/subsets/test/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    output:
        html_report = f"work_folder{pn}/subsets/report/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}.html"
    shell:
        """
        Rscript -e "rmarkdown::render(
          '{params.rmd_location}',
          output_file = '{output.html_report}',
          params = list(
            full_pos_file = '{input.full_pos}',
            full_neg_file = '{input.full_neg}',
            train_pos_file = '{input.train_pos}',
            train_neg_file = '{input.train_neg}',
            train_maxflow_pos_file = '{input.train_maxflow_pos}',
            train_maxflow_neg_file = '{input.train_maxflow_neg}',
            validate_pos_file = '{input.validate_pos}',
            validate_neg_file = '{input.validate_neg}',
            validate_maxflow_pos_file = '{input.validate_maxflow_pos}',
            validate_maxflow_neg_file = '{input.validate_maxflow_neg}',
            test_pos_file = '{input.test_pos}',
            test_neg_file = '{input.test_neg}',
            test_maxflow_pos_file = '{input.test_maxflow_pos}',
            test_maxflow_neg_file = '{input.test_maxflow_neg}'
          )
        )"
        """