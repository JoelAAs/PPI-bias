import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import pymc as pm
import ray


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

        trace = pm.sample(samples, tune=tunings, chains=4, cores=1, return_inferencedata=True, progressbar=True)

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).values

    ordered_values = [val for pair in zip(beta_detection_mu, beta_detection_sd) for val in pair]

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    return ordered_values + [beta_bait_mu, beta_bait_sd, n_diverging]

@ray.remote
def process_prey(prey_df, cl_categories, bait_name, prey_name, detection_matrix):
    bait_matrix = np.zeros((prey_df["n_tested"].sum(), len(cl_categories) + 2))
    bait_matrix[:,1] = 1
    c_row = 0
    for i, row in prey_df.iterrows():
        interaction_observations = row['n_observed']
        n_studies_bait = row['n_tested']
        cl_id = row["cl_id"]
        cl_match = np.where(cl_categories == cl_id)[0] ## try except key error
        bait_matrix[c_row:(c_row+n_studies_bait), cl_match+2] = 1
        bait_matrix[c_row:(c_row+interaction_observations), 0] = 1
        c_row += n_studies_bait

    interaction_observations_total = prey_df['n_observed'].sum()
    n_studies_bait_total = prey_df['n_tested'].sum()
    input_matrix = np.concat((detection_matrix, bait_matrix))

    detection_parameters = detection_model(input_matrix)

    # Return a tab-separated string for writing later
    return f"{bait_name}\t{prey_name}\t{interaction_observations_total}\t{n_studies_bait_total}\t" + \
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
    cl_categories = np.array(["CVCL_0030", "CVCL_0063", "CVCL_0291"])

    def bin_it(list, size):
        if size > len(list):
            return [list]
        return [list[:size]] + bin_it(list[size:], size)

    batched_prey = bin_it(tested_prey, 200)
    for preys in batched_prey:
        start = datetime.now()

        ray.init(num_cpus=args.workers)
        futures = []
        for prey in preys:
            baseline_pod_prey     = baseline_pod[baseline_pod["gene_name"] == prey]
            detection_matrix      = get_one_hot_untargeted(cl_categories, baseline_pod_prey["cl_id"])
            detection_matrix[:,0] =  baseline_pod_prey["value"]

            ppi_ss = bait_df[bait_df["gene_name_prey"] == prey]
            bait_name = ppi_ss["gene_name_bait"].unique()[0] # allways the same, do better
            futures.append(process_prey.remote(
                prey_df=ppi_ss,
                cl_categories=cl_categories,
                bait_name=bait_name,
                prey_name=prey,
                detection_matrix=detection_matrix))

        results = ray.get(futures)

        stop = datetime.now()
        print(f"Estimated probability for {len(preys)} in {stop - start} acorss {args.workers} cores")

        with open(args.bait_output, "a+") as w:
            for line in results:
                w.write(line + "\n")

        ray.shutdown()

if __name__ == "__main__":
    main()