import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import pymc as pm
import ray


def get_one_hot_untargeted(cl_categories, cell_line_annotations):
    untargeted_matrix = np.zeros(shape=(len(cell_line_annotations), len(cl_categories) + 2))
    for i, cl_annotation in enumerate(cell_line_annotations):
        cl_match = np.where(cl_categories == cl_annotation)[0]  ## try except key error

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
        # start = pm.find_MAP()
        beta_detection = pm.Normal('beta_detection', mu=0, sigma=10,
                                   shape=value_matrix.shape[1] - 2)
        beta_bait = pm.Normal('beta_bait', mu=0, sigma=10)

        logit_p = pm.math.dot(value_matrix[:, 2:], beta_detection) + beta_bait * value_matrix[:, 1]
        p = pm.Deterministic('p', pm.math.sigmoid(logit_p))

        y_obs = pm.Bernoulli('y_obs', p=p, observed=value_matrix[:, 0])

        trace = pm.sample(samples, tune=tunings, chains=4, cores=1, target_accept=.9, return_inferencedata=True,
                          progressbar=False)

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).values

    ordered_values = [val for pair in zip(beta_detection_mu, beta_detection_sd) for val in pair]

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    return ordered_values + [beta_bait_mu, beta_bait_sd, n_diverging]


@ray.remote
def process_prey(interaction_row, obs_c, tested_c, cl_categories, detection_matrix):
    bait_matrix = np.zeros((interaction_row[tested_c].sum(), len(cl_categories) + 2))
    bait_matrix[:, 1] = 1
    i = 0
    cl_i = 2
    for cl_obs, cl_test in zip(obs_c, tested_c):
        bait_matrix[i:(i + interaction_row[cl_obs]), 0] = 1
        bait_matrix[i:(i + interaction_row[cl_test]), cl_i] = 1
        i += interaction_row[cl_test]
        cl_i += 1
    input_matrix = np.concat((detection_matrix, bait_matrix))

    detection_parameters = detection_model(input_matrix)

    # Return a tab-separated string for writing later
    return "\t".join(map(str, interaction_row.values)) + "\t" + "\t".join(map(str, detection_parameters))


def main():
    parser = argparse.ArgumentParser(description="Process bait interaction data with multiprocessing.")
    parser.add_argument("--prey_tested", required=True, help="Path to prey combination TSV file")
    parser.add_argument("--pod_base_reform", required=True, help="Path to baseline pod TSV file")
    parser.add_argument("--bait_output", required=True, help="Path to output file for bait parameters")
    parser.add_argument("--workers", type=int,
                        help="Number of worker processes for ray")
    args = parser.parse_args()
    prey_interaction_df = pd.read_csv(args.prey_tested, sep="\t")
    numeric_cols = prey_interaction_df.columns[1:]
    prey_interaction_df[numeric_cols] = prey_interaction_df[numeric_cols].astype(int)
    baseline_pod = pd.read_csv(args.pod_base_reform, sep="\t")

    cl_categories = np.array(["CVCL_0030", "CVCL_0063", "CVCL_0291"])

    batch_size = 200
    i = 0
    while prey_interaction_df.shape[0] > i:
        def get_order(cols, pattern, pos=0):
            if pattern in cols[0]:
                return pos
            else:
                return get_order(cols[1:], pattern, pos + 1)

        observed_cols = np.array([col_name for col_name in prey_interaction_df.columns if "n_observed" in col_name])
        idx_o = [get_order(observed_cols, cl_order) for cl_order in cl_categories]
        observed_cols = observed_cols[idx_o]

        tested_cols = np.array([col_name for col_name in prey_interaction_df.columns if "n_tested" in col_name])
        idx_t = [get_order(tested_cols, cl_order) for cl_order in cl_categories]
        tested_cols = tested_cols[idx_t]

        start = datetime.now()

        ray.init(num_cpus=args.workers)
        futures = []
        for _, row in prey_interaction_df.iloc[i:(i + batch_size)].iterrows():
            prey = row["gene_name_prey"]
            baseline_pod_prey = baseline_pod[baseline_pod["gene_name"] == prey]
            detection_matrix = get_one_hot_untargeted(cl_categories, baseline_pod_prey["cl_id"])
            detection_matrix[:, 0] = baseline_pod_prey["value"]

            futures.append(process_prey.remote(
                interaction_row=row,
                obs_c=observed_cols,
                tested_c=tested_cols,
                cl_categories=cl_categories,
                detection_matrix=detection_matrix))

        results = ray.get(futures)

        stop = datetime.now()
        print(f"Estimated probability for {batch_size} in {stop - start} acorss {args.workers} cores")

        with open(args.bait_output, "a+") as w:
            for line in results:
                w.write(line + "\n")

        ray.shutdown()


if __name__ == "__main__":
    main()
