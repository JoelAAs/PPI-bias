
def get_expected_input(wc):
    if wc.dataset == "goldensplit":
        data = f"data_{wc.network_type}"
        selection = wc.dataset
    else:
        if wc.network_type == "directional":
            selection = "maxflow"
        elif wc.network_type == "undirectional":
            selection="undirectionalbalanced"
        else:
            raise ValueError(f"unknown network type {wc.network_type}")

        pos_limit = config["models"][wc.model_configuration]["pos"]
        neg_limit = config["models"][wc.model_configuration]["neg"]

        if re.search("-random",wc.partition):
            pos_limit = config["models"][wc.model_configuration]["pos"]
            posdata =  f"{wc.dataset}_{wc.network_type}_limit_{neg_limit}_poslim_{pos_limit}_{wc.partition.split("-")[0]}"
            negdata =  f"{wc.dataset}_{wc.network_type}_limit_{neg_limit}_poslim_{pos_limit}_{wc.partition}"
            return [
                f"work_folder{pn}/subsets/train/{selection}/{posdata}_pos.pq",
                f"work_folder{pn}/subsets/train/{selection}/{negdata}_neg.pq",
                f"work_folder{pn}/subsets/validation/{selection}/{posdata}_pos.pq",
                f"work_folder{pn}/subsets/validation/{selection}/{negdata}_neg.pq",
                f"work_folder{pn}/subsets/test/{selection}/{posdata}_pos.pq",
                f"work_folder{pn}/subsets/test/{selection}/{negdata}_neg.pq"
            ]

        else:
            data =  f"{wc.dataset}_{wc.network_type}_limit_{neg_limit}_poslim_{pos_limit}_{wc.partition}"

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
        params =      f"work_folder{pn}/classification/randomforest/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_model_parameters.txt",
        saved_model = f"work_folder{pn}/classification/randomforest/model/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}_model_parameters.joblib"
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
            --saved_model {output.saved_model} 
        """



rule get_tvt_degree_balance:
    params:
        script_location = "PPIClassification/ModelEvaluation/degree_balance_metrics.py"
    input:    
        data = lambda wc: get_expected_input(wc)
    output:
        degree_balance = f"work_folder{pn}/subsets/degree_balance/{{dataset}}_{{network_type}}_{{model_configuration}}_{{partition}}.csv"
    shell:
        """
        python3 {params.script_location} \
            --pos_train {input.data[0]} \
            --neg_train {input.data[1]} \
            --pos_val {input.data[2]} \
            --neg_val {input.data[3]} \
            --pos_test {input.data[4]} \
            --neg_test  {input.data[5]} \
            --output_file {output.degree_balance} \
            --network_type {wildcards.network_type}
        """
