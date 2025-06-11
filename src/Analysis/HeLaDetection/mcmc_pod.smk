import pandas as pd
import glob
import os


def get_bait_parameters(wildcards):
    BAIT_FOLDER = checkpoints.estimate_bait_interaction.get().output[0]
    bait_files = glob.glob(BAIT_FOLDER +"/*")

    expected = [
        c_file.replace(".csv", "_parameters.csv") for c_file in bait_files
    ]
    return expected



rule split_double_columns:
    input:
        pod_base = "data/absent_0_present_1_selected_N07444_M04547.csv"
    output:
        pod_base_reform = "work_folder/POD/HeLa.csv"
    run:
        baseline_pod = pd.read_csv(input.pod_base).fillna(0)
        cols_to_split = [col for col in baseline_pod.columns if ';' in col]

        for col in cols_to_split:
            new_cols = col.split(';')
            for new_col in new_cols:
                baseline_pod[new_col] = baseline_pod[col] # TODO: DataFrame is highly fragmented ... fix
            baseline_pod = baseline_pod.drop(columns=[col])

        baseline_pod = baseline_pod.drop(columns=["Unnamed: 0"]) # R index prob
        baseline_pod.to_csv(output.pod_base_reform, sep="\t", index=False)


checkpoint estimate_bait_interaction:
    params:
        min_tested = 4,
        cellines = ["CVCL_0030", "CVCL_0291", "CVCL_0063"]
    input:
        pod_base_reform = "~/resources/proteome/cl_proteome/shared_detection.csv",
        cl_specific_interactions = "data/CL_annotated_bait_prey.csv"
    output:
        bait_folder = directory("work_folder/analysis/Hela_pod/baits")
    run:
        os.makedirs(output.bait_folder, exist_ok=True)
        PPI_interactions = pd.read_csv(input.cl_specific_interactions, sep="\t")
        baseline_pod = pd.read_csv(input.pod_base_reform, sep="\t")

        PPI_interactions = PPI_interactions[
            (PPI_interactions["gene_name_bait"] != PPI_interactions["gene_name_prey"]) &
            (PPI_interactions["cl_id"].isin(params.cellines))
        ]
        PPI_interactions["study_id"] = PPI_interactions[
            ["pubmed_id", "detection_method"]
        ].apply(lambda row: "-".join(map(str, row)), axis=1)
        n_studies = PPI_interactions.groupby(["gene_name_bait", "cl_id"], as_index=False
        )["study_id"].nunique()
        n_studies = n_studies.rename(
            {
                "study_id": "n_tested"
            }, axis = 1
        )
        n_bait_prey_tests = PPI_interactions.groupby(
            ["gene_name_bait", "gene_name_prey", "cl_id"])["study_id"].nunique()

        total_studies = n_studies.groupby("gene_name_bait")["n_tested"].sum()
        baits_to_keep = total_studies[total_studies > params.min_tested].index.values
        n_studies = n_studies[n_studies["gene_name_bait"].isin(baits_to_keep)]

        for bait_name in n_studies["gene_name_bait"].unique():
            cl_ss = n_studies[n_studies["gene_name_bait"] == bait_name]
            with open(f"{output.bait_folder}/{bait_name}.csv","w") as w:
                w.write("\t".join(
                    [
                        "gene_name_bait",
                        "gene_name_prey",
                        "n_observed",
                        "n_tested",
                        "cl_id"
                    ]
                ) + "\n")
                for i, (_, cl_id, n_studies_bait) in cl_ss.iterrows():
                    for possible_prey in baseline_pod["gene_name"].unique():
                        try:
                            interaction_observations = n_bait_prey_tests[bait_name][possible_prey][cl_id]
                        except KeyError:
                            interaction_observations = 0

                        w.write("\t".join(
                            map(str,[
                                bait_name,
                                possible_prey,
                                interaction_observations,
                                n_studies_bait,
                                cl_id
                            ])
                        ) + "\n")


rule fit_parameters:
        params:
            workers = 25
        input:
            bait = "work_folder/analysis/Hela_pod/baits/{bait}.csv",
            pod_base_reform = "~/resources/proteome/cl_proteome/shared_detection.csv"
        output:
            bait_parameters = "work_folder/analysis/Hela_pod/baits/{bait}_parameters.csv"
        shell:
            """
            python src/Analysis/HeLaDetection/fit_model_mp.py \
                --bait {input.bait} \
                --pod_base_reform {input.pod_base_reform} \
                --bait_output {output.bait_parameters} \
                --workers {params.workers}
            """



rule aggregate:
        input:
            bait_parameters = get_bait_parameters
        output:
            aggregate_parameters = "work_folder/analysis/Hela_pod/all_parameters.csv"
        run:
            with open(output.aggregate_parameters, "w") as w:
                w.write("\t".join(
                    [
                        "gene_name_bait",
                        "gene_name_prey",
                        "n_observed",
                        "n_tested",
                        "beta_prediction_mean",
                        "beta_prediction_sd",
                        "beta_bait_mean",
                        "beta_bait_sd",
                        "n_divergences"
                    ]
                ) + "\n")
                for bait_parameters in input.bait_parameters:
                    with open(bait_parameters, "r") as f:
                        for line in f:
                            w.write(line)


