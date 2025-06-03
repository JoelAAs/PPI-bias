import argparse
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os
import pymc as pm

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

def process_row(row, baseline_pod):
    bait_name = row['gene_name_prey']
    possible_prey = row['gene_name_prey']
    interaction_observations = row['n_observed']
    n_studies_bait = row['n_tested']

    detections = baseline_pod[possible_prey].tolist()
    baits = [0] * len(detections)

    detections += [1] * interaction_observations + [0] * (n_studies_bait - interaction_observations)
    baits += [1] * n_studies_bait

    detection_parameters = detection_model(detections, baits)

    # Return a tab-separated string for writing later
    return f"{bait_name}\t{possible_prey}\t{interaction_observations}\t{n_studies_bait}\t" + \
           "\t".join(map(str, detection_parameters))

def main():
    parser = argparse.ArgumentParser(description="Process bait interaction data with multiprocessing.")
    parser.add_argument("--bait", required=True, help="Path to bait input TSV file")
    parser.add_argument("--pod_base_reform", required=True, help="Path to baseline pod TSV file")
    parser.add_argument("--bait_output", required=True, help="Path to output file for bait parameters")
    parser.add_argument("--workers", type=int,
                        help="Number of worker processes for multiprocessing")

    args = parser.parse_args()

    bait_df = pd.read_csv(args.bait, sep="\t")
    baseline_pod = pd.read_csv(args.pod_base_reform, sep="\t")
    rows = [row for _, row in bait_df.iterrows()]

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        from functools import partial
        func = partial(process_row, baseline_pod=baseline_pod)
        results = list(executor.map(func, rows))

    with open(args.bait_output, "w") as w:
        for line in results:
            w.write(line + "\n")

if __name__ == "__main__":
    main()