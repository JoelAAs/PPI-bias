import concurrent.futures
import pandas as pd
import pymc as pm
import glob
import os


def model_per_bait(output_csv, baseline_pod, n_bait_prey_tests, n_studies_bait, bait_name):
    with open(output_csv,"w") as w:
        for possible_prey in baseline_pod.columns.values:

            detections = baseline_pod[possible_prey].tolist()
            baits = [0] * len(detections)

            try:
                interaction_observations = n_bait_prey_tests[bait_name][possible_prey]
            except KeyError:
                interaction_observations = 0

            detections += [1] * interaction_observations + [0] * (n_studies_bait - interaction_observations)
            baits += [1] * n_studies_bait

            detection_parameters = detection_model(detections,baits)

            w.write(f"{bait_name}\t{possible_prey}"
                    f"\t{interaction_observations}"
                    f"\t{n_studies_bait}" + "\t".join(map(str,detection_parameters)) + "\n")

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


rule estimate_bait_interaction:
    params:
        prior = 1,
        min_tested = 4,
        celline = "CVCL_0030",
        max_workers = 4
    input:
        pod_base = "data/absent_0_present_1_selected_N07444_M04547.csv",
        cl_specific_interactions = "data/CL_annotated_bait_prey.csv"
    output:
        bait_folder = directory("work_folder/analysis/Hela_pod/baits"),
        all_found   = "work_folder/analysis/Hela_pod/logistic_mcmc.csv"
    run:
        os.makedirs(output.bait_folder, exist_ok=True)
        PPI_interactions = pd.read_csv(input.cl_specific_interactions, sep="\t")
        baseline_pod = pd.read_csv(input.pod_base).fillna(0)
        cols_to_split = [col for col in baseline_pod.columns if ';' in col]

        for col in cols_to_split:
            new_cols = col.split(';')
            for new_col in new_cols:
                baseline_pod[new_col] = baseline_pod[col] # TODO: DataFrame is highly fragmented ... fix
            baseline_pod = baseline_pod.drop(columns=[col])

        baseline_pod= baseline_pod.drop(columns=["Unnamed: 0"]) # R index prob

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

        with concurrent.futures.ThreadPoolExecutor(max_workers=params.max_workers) as executor:
            futures = [
                executor.submit(model_per_bait,
                    f"{output.bait_folder}/{bait_name}.csv",
                    baseline_pod, n_bait_prey_tests,
                    n_studies_bait, bait_name)
                for _, (bait_name, n_studies_bait) in n_studies.iterrows()
            ]

            for future in concurrent.futures.as_completed(futures):
                future.result()

        with open(output.all_found, "w") as w:
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
            for (bait_name, n_studies_bait) in n_studies.iterrows():
                with open(f"{output.bait_folder}/{bait_name}.csv", "r") as f:
                    for line in f:
                        w.write(line)