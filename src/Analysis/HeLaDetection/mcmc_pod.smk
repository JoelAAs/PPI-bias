import pandas as pd
import pymc as pm
import glob
import os


def get_bait_parameters(wildcards):
    BAIT_FOLDER = checkpoints.estimate_bait_interaction.get().output[0]
    bait_files = glob.glob(BAIT_FOLDER +"/*")

    expected = [
        c_file.replace(".csv", "_parameters.csv") for c_file in bait_files
    ]
    return expected

def detection_model(y_detections, x_baits, samples=1000, tunings=500, cores=2):
    with pm.Model() as logistic_model:
        beta_detection = pm.Normal('beta_detection',mu=0,sigma=10)
        beta_bait = pm.Normal('beta_bait',mu=0,sigma=10)


        logit_p = beta_detection + beta_bait * x_baits

        p = pm.Deterministic('p',pm.math.sigmoid(logit_p))
        y_obs = pm.Bernoulli('y_obs',p=p,observed=y_detections)

        trace = pm.sample(samples,tune=tunings,cores=cores,return_inferencedata=True)

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).item()
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).item()

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()

    return beta_detection_mu, beta_detection_sd, beta_bait_mu, beta_bait_sd


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
        prior = 1,
        min_tested = 4,
        celline = "CVCL_0030",
        max_workers = 4
    input:
        pod_base_reform = "work_folder/POD/HeLa.csv",
        cl_specific_interactions = "data/CL_annotated_bait_prey.csv"
    output:
        bait_folder = directory("work_folder/analysis/Hela_pod/baits")
    run:
        os.makedirs(output.bait_folder, exist_ok=True)
        PPI_interactions = pd.read_csv(input.cl_specific_interactions, sep="\t")
        baseline_pod = pd.read_csv(input.pod_base_reform, sep="\t")

        PPI_interactions = PPI_interactions[
            (PPI_interactions["gene_name_bait"] != PPI_interactions["gene_name_prey"]) &
            (PPI_interactions["cl_id"] == params.celline)
        ]
        PPI_interactions["study_id"] = PPI_interactions[
            ["pubmed_id", "detection_method"]
        ].apply(lambda row: "-".join(map(str, row)), axis=1)
        n_studies = PPI_interactions.groupby("gene_name_bait", as_index=False
        )["study_id"].nunique()
        n_studies = n_studies.rename(
            {
                "study_id": "n_tested"
            }, axis = 1
        )
        n_bait_prey_tests = PPI_interactions.groupby(
            ["gene_name_bait", "gene_name_prey"])["study_id"].nunique()

        n_studies = n_studies[n_studies["n_tested"] >= params.min_tested]

        for i, (bait_name, n_studies_bait) in n_studies.iterrows():
            with open(f"{output.bait_folder}/{bait_name}.csv","w") as w:
                for possible_prey in baseline_pod.columns.values:
                    try:
                        interaction_observations = n_bait_prey_tests[bait_name][possible_prey]
                    except KeyError:
                        interaction_observations = 0

                    w.write("\t".join(
                        map(str,[
                            bait_name,
                            possible_prey,
                            interaction_observations,
                            n_studies_bait
                        ])
                    ) + "\n")


rule fit_parameters:
        input:
            bait = "work_folder/analysis/Hela_pod/baits/{bait}.csv",
            pod_base_reform = "work_folder/POD/HeLa.csv"
        output:
            bait_parameters = "work_folder/analysis/Hela_pod/baits/{bait}_parameters.csv"
        run:
            bait_df = pd.read_csv(input.bait, sep="\t", header=False)
            baseline_pod = pd.read_csv(input.pod_base_reform, sep="\t")
            with open(output.bait_parameters, "w") as w:
                for i, (bait_name, possible_prey, interaction_observations, n_studies_bait) in bait_df.iterrows():

                    detections = baseline_pod[possible_prey].tolist()
                    baits = [0] * len(detections)

                    detections += [1] * interaction_observations + [0] * (
                            n_studies_bait - interaction_observations)
                    baits += [1] * n_studies_bait

                    detection_parameters = detection_model(detections,baits)

                    w.write(f"{bait_name}\t{possible_prey}"
                            f"\t{interaction_observations}"
                            f"\t{n_studies_bait}" + "\t".join(map(str,detection_parameters)) + "\n")


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
                    ]
                ) + "\n")
                for bait_parameters in input.bait_parameters:
                    with open(bait_parameters, "r") as f:
                        for line in f:
                            w.write(line)


