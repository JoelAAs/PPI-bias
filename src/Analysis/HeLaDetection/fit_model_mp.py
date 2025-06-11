import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import pymc as pm
import ray

# import logging
# logger = logging.getLogger("pymc")
# logger.propagate = False

def get_input_matrix(cl_categories):



def get_one_hot_untargeted(cl_categories, cell_line_annotations):
    untargeted_matrix = np.zeros(shape = (len(cell_line_annotations), len(cl_categories) + 2))
    for i, cl_annotation in enumerate(cell_line_annotations):
        cl_match = np.where(cl_categories == cl_annotation)[0] ## try except key error

        untargeted_matrix[i, cl_match + 2] = 1

    return untargeted_matrix

def detection_model(value_matrix, samples=1000, tunings=500):
    """
    :param value_matrix: shape (detection, bait, cell_line variables)
    :param samples: number of MCMC sampling
    :param tunings: number of burn-in samples
    :return:
    """
    with pm.Model() as multi_env_model:
        beta_detection = pm.Normal('beta_detection', mu=0, sigma=10,
                                   shape=value_matrix.shape[1] - 2)
        beta_bait = pm.Normal('beta_bait', mu=0, sigma=10)

        logit_p = pm.math.dot(value_matrix[:,2:], beta_detection) + beta_bait * value_matrix[:,1]
        p = pm.Deterministic('p', pm.math.sigmoid(logit_p))

        y_obs = pm.Bernoulli('y_obs', p=p, observed=value_matrix[:,0])

        trace = pm.sample(samples, tune=tunings, cores=1, return_inferencedata=True, progressbar=False)

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).item()
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).item()

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    return beta_detection_mu, beta_detection_sd, beta_bait_mu, beta_bait_sd, n_diverging

@ray.remote
def process_row(prey_df, prey_baseline_pod):

    
    bait_name = row['gene_name_bait']
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

    tested_prey = bait_df["gene_name_prey"].unique()



    def bin_it(list, size):
        if size > len(list):
            return [list]
        return [list[:size]] + bin_it(list[size:], size)

    batched_rows = bin_it(all_rows, 200)
    for rows in batched_rows:
        start = datetime.now()

        ray.init(num_cpus=args.workers)

        futures = [process_row.remote(row, baseline_pod.copy()) for row in rows]

        results = ray.get(futures)

        stop = datetime.now()
        print(f"Estimated probability for {len(rows)} in {stop - start} acorss {args.workers} cores")

        with open(args.bait_output, "a+") as w:
            for line in results:
                w.write(line + "\n")

        ray.shutdown()

if __name__ == "__main__":
    main()