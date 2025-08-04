import argparse
import time
from datetime import datetime
import os
import random
os.environ["RAY_DEDUP_LOGS"] = "0"
import numpy as np
import pandas as pd
import ray

def get_one_hot_untargeted(cl_categories, untargeted_protein_relactive_abundance):
    untargeted_matrix = np.zeros(shape=(len(untargeted_protein_relactive_abundance), len(cl_categories) + 2))
    i=0
    for _, (abundance, cl_annotation) in untargeted_protein_relactive_abundance.iterrows():
        cl_match = np.where(cl_categories == cl_annotation)[0]  ## try except key error
        untargeted_matrix[i, cl_match + 2] = 1
        untargeted_matrix[i, 0] = int(abundance != 0)
        i += 1
    return untargeted_matrix


def get_prey_sd_mu(df, prey):
    df_mu = df.groupby("cell_line", as_index=False).mean().rename({prey: "mu"}, axis=1)
    df_sd = df.groupby("cell_line", as_index=False).std().rename({prey:"sd"}, axis=1)
    return df_mu.merge(df_sd, on="cell_line").set_index("cell_line")


@ray.remote(num_cpus=1)
def process_prey_pod(interaction_row, obs_c, tested_c, cl_categories, detection_matrix, samples=1000, tunings=500):
    import uuid
    import pymc as pm
    import logging
    logger = logging.getLogger("pymc")
    logger.setLevel(logging.ERROR)

    start = datetime.now()

    N_tests = interaction_row[tested_c].sum()
    bait_matrix = np.zeros((N_tests, len(cl_categories) + 2))
    bait_matrix[:, 1] = 1
    i = 0
    cl_i = 2

    for cl_obs, cl_test in zip(obs_c, tested_c):
        bait_matrix[i:(i + interaction_row[cl_obs]), 0] = 1
        bait_matrix[i:(i + interaction_row[cl_test]), cl_i] = 1
        i += interaction_row[cl_test]
        cl_i += 1

    value_matrix = np.concat((detection_matrix, bait_matrix))

    with pm.Model() as multi_env_model:

        beta_detection = pm.Normal('beta_detection', mu=0, sigma=10, shape=value_matrix.shape[1] - 2)
        beta_bait = pm.Normal('beta_bait', mu=0, sigma=10)

        logit_p = pm.math.dot(value_matrix[:, 2:], beta_detection) + beta_bait * value_matrix[:, 1]
        p = pm.Deterministic('p', pm.math.sigmoid(logit_p))

        y_obs = pm.Bernoulli('y_obs', p=p, observed=value_matrix[:, 0])


        try:
            trace = pm.sample(samples, tune=tunings, chains=4, cores=1, target_accept=0.95, return_inferencedata=True,
                              progressbar=False)
        except (TimeoutError, AssertionError):
            print("failed, redoing")
            time.sleep(15) # If multiple jobs tries to compile. Doesn't solve it just makes it less likely
            trace = pm.sample(samples, tune=tunings, chains=4, cores=1, target_accept=0.95, return_inferencedata=True,
                              progressbar=False)

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).values

    ordered_values = [val for pair in zip(beta_detection_mu, beta_detection_sd) for val in pair]

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    detection_parameters =  ordered_values + [beta_bait_mu, beta_bait_sd, n_diverging]
    end = datetime.now()
    print(f"A job took {end-start} time")
    # Return a tab-separated string for writing later
    return "\t".join(map(str, interaction_row.values)) + "\t" + "\t".join(map(str, detection_parameters))



@ray.remote(num_cpus=1)
def process_prey_abundance(interaction_row, obs_c, tested_c, cl_categories, prey_distribution, samples=1500, tunings=1000):
    import pymc as pm
    import logging
    logger = logging.getLogger("pymc")
    logger.setLevel(logging.ERROR)

    start = datetime.now()
    target_accept = 0.8
    N_tests = interaction_row[tested_c].sum()
    bait_matrix = np.zeros((N_tests, len(cl_categories) + 2))
    bait_matrix[:, 1] = 1
    i = 0
    cl_i = 2
    cl_idx = np.zeros(bait_matrix.shape[0], dtype=int)

    for cl_obs, cl_test in zip(obs_c, tested_c):
        cl_idx[i:(i + interaction_row[cl_test])] = cl_i-2
        bait_matrix[i:(i + interaction_row[cl_obs]), 0] = 1
        bait_matrix[i:(i + interaction_row[cl_test]), cl_i] = 1
        i += interaction_row[cl_test]
        cl_i += 1

    ## lets do normal distribution so, therefore do log2
    mu_untargeted = np.zeros(len(cl_categories))
    sd_untargeted = np.zeros(len(cl_categories))
    for cl, (mu, sd) in prey_distribution.iterrows():
        mu_untargeted[np.where(cl_categories == cl)[0]] = mu
        sd_untargeted[np.where(cl_categories == cl)[0]] = sd

    target_accept_low = True
    non_div_run = False
    while not non_div_run:
        with pm.Model() as multi_env_model:
            b_bait = pm.Normal('b_bait', mu=0, sigma=10)
            x_untargeted = pm.Normal(
                'x_untargeted',
                mu=mu_untargeted,
                sigma=sd_untargeted,
                shape=mu_untargeted.shape[0]
            )

            x_to_samples = x_untargeted[cl_idx]
            p_y = pm.Deterministic("mu_y", pm.math.sigmoid(b_bait*(2**x_to_samples))) # we estimate mean on log
            y_lh = pm.Bernoulli('y_lh', p=p_y, observed=bait_matrix[:, 0])

            try:
                trace = pm.sample(samples, tune=tunings, chains=4, cores=1, target_accept=target_accept, return_inferencedata=True,
                                  progressbar=False)
            except (TimeoutError, AssertionError):
                print("failed, redoing")
                time.sleep(15) # If multiple jobs tries to compile. Doesn't solve it just makes it less likely
                trace = pm.sample(samples, tune=tunings, chains=4, cores=1, target_accept=target_accept, return_inferencedata=True,
                                  progressbar=False)
        n_diverging = sum(sum(trace.sample_stats.diverging.values))
        if not target_accept_low:
            non_div_run = True
        elif n_diverging > 0:
            target_accept = .99
            target_accept_low = False
        else:
            non_div_run = True
            

    beta_detection_mu = trace.posterior["x_untargeted"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["x_untargeted"].std(("chain", "draw")).values

    ordered_values = [val for pair in zip(beta_detection_mu, beta_detection_sd) for val in pair]

    beta_bait_mu = trace.posterior["b_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["b_bait"].std(("chain", "draw")).item()

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    detection_parameters =  ordered_values + [beta_bait_mu, beta_bait_sd, n_diverging]
    end = datetime.now()
    print(f"A job took {end-start} time at target_accept {target_accept}")
    # Return a tab-separated string for writing later
    return "\t".join(map(str, interaction_row.values)) + "\t" + "\t".join(map(str, detection_parameters))


def main():
    parser = argparse.ArgumentParser(description="Process bait interaction data with multiprocessing.")
    parser.add_argument("--prey_tested", required=True, help="Path to prey combination TSV file")
    parser.add_argument("--abundance_cell_lines", required=True, help="Path to baseline pod TSV file")
    parser.add_argument("--bait_output", required=True, help="Path to output file for bait parameters")
    parser.add_argument("--abundance", type=int, help="Set to 1 for abundance, 0 for probability of detection")
    parser.add_argument("--workers", type=int, help="Number of worker processes for ray")
    parser.add_argument("--batch_size", type=int, help="Number of tests performed before writing to file")
    args = parser.parse_args()
    prey_interaction_df = pd.read_csv(args.prey_tested, sep="\t")
    numeric_cols = prey_interaction_df.columns[1:]
    prey_interaction_df[numeric_cols] = prey_interaction_df[numeric_cols].astype(int)
    cell_line_abundance = pd.read_csv(args.abundance_cell_lines, sep="\t")
    cl_categories = np.array(["CVCL_0030", "CVCL_0063", "CVCL_0291"])

    batch_size = args.batch_size
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
        for j, row in prey_interaction_df.iloc[i:(i + batch_size)].iterrows():
            ray_task_env = {"PYTENSOR_FLAGS": f"compiledir=/tmp/compiledir_{j}"}
            prey = row["gene_name_prey"]
            if args.abundance == 1:
                prey_mu_sd = get_prey_sd_mu(cell_line_abundance[[prey, "cell_line"]].dropna(), prey)
                futures.append(process_prey_abundance.options(
                    runtime_env={"env_vars": ray_task_env}
                ).remote(
                    interaction_row=row,
                    obs_c=observed_cols,
                    tested_c=tested_cols,
                    cl_categories=cl_categories,
                    prey_distribution=prey_mu_sd))

            elif args.abundance == 0:
                detection_matrix = get_one_hot_untargeted(
                   cl_categories,
                   cell_line_abundance[[prey, "cell_line"]].fillna(0)) #remove non-observations

                futures.append(process_prey_pod.options(
                    runtime_env={"env_vars": ray_task_env}
                ).remote(
                    interaction_row=row,
                    obs_c=observed_cols,
                    tested_c=tested_cols,
                    cl_categories=cl_categories,
                    detection_matrix=detection_matrix))

            elif args.abundance:
                raise ValueError(f"Unknown input value {args.abundance} for abundance")

        results = ray.get(futures)
        stop = datetime.now()
        print(f"Estimated probability for {batch_size} in {stop - start} acorss {args.workers} cores")

        with open(args.bait_output, "a+") as w:
            for line in results:
                w.write(line + "\n")

        ray.shutdown()
        i += batch_size


if __name__ == "__main__":
    main()
