import argparse
import time
from datetime import datetime
import os
import arviz as az
import math
os.environ["RAY_DEDUP_LOGS"] = "0"
import numpy as np
import pandas as pd
import ray


def get_one_hot_untargeted(cl_categories, untargeted_protein_relactive_abundance):
    untargeted_matrix = np.zeros(shape=(len(untargeted_protein_relactive_abundance), len(cl_categories) + 2))
    i = 0
    for _, (abundance, cl_annotation) in untargeted_protein_relactive_abundance.iterrows():
        cl_match = np.where(cl_categories == cl_annotation)[0]  ## try except key error
        untargeted_matrix[i, cl_match + 2] = 1
        untargeted_matrix[i, 0] = int(abundance != 0)
        i += 1
    return untargeted_matrix


def get_prey_sd_mu(df, prey):
    df_mu = df.groupby("cell_line", as_index=False).mean().rename({prey: "mu"}, axis=1)
    df_sd = df.groupby("cell_line", as_index=False).std().rename({prey: "sd"}, axis=1)
    return df_mu.merge(df_sd, on="cell_line").set_index("cell_line")


@ray.remote(num_cpus=1)
def process_prey_pod(interaction_row, obs_c, tested_c, cl_categories, detection_matrix, j, samples=1000, tunings=500):
    import pymc as pm
    import logging
    from pymc.variational.callbacks import CheckParametersConvergence

    logger = logging.getLogger("pymc")
    logger.setLevel(logging.ERROR)

    start = datetime.now()
    target_accept = .8
    n_tests = interaction_row[tested_c].sum()
    bait_matrix = np.zeros((n_tests, len(cl_categories) + 2))
    bait_matrix[:, 1] = 1
    i = 0
    cl_i = 2

    for cl_obs, cl_test in zip(obs_c, tested_c):
        bait_matrix[i:(i + interaction_row[cl_obs]), 0] = 1
        bait_matrix[i:(i + interaction_row[cl_test]), cl_i] = 1
        i += interaction_row[cl_test]
        cl_i += 1

    value_matrix = np.concat((detection_matrix, bait_matrix))
    target_accept_low = True
    non_div_run = False
    n_chains = 4

    while not non_div_run:
        with pm.Model() as multi_env_model:
            logit = lambda x: -math.log(1 / x - 1)
            sd_prior = (logit(.99) - logit(.25)) / (1.96 * 2)  # due to filtering the lowest pod would be 75 %
            mu_prior = logit(0.62)
            beta_detection = pm.Normal('beta_detection', mu=mu_prior, sigma=sd_prior, shape=value_matrix.shape[1] - 2)
            beta_bait = pm.SkewNormal('beta_bait', mu=0, sigma=10, alpha=0)

            logit_p = pm.math.dot(value_matrix[:, 2:], beta_detection) + beta_bait * value_matrix[:, 1]

            lp = pm.Deterministic('p', logit_p)
            y_obs = pm.Bernoulli('y_obs', logit_p=lp, observed=value_matrix[:, 0])
            vi_approx = pm.fit(
                method="advi",
                callbacks=[CheckParametersConvergence(diff="absolute")])
            init_values = []
            for _ in range(n_chains):
                vi_samples = vi_approx.sample(1)
                for var_name in vi_samples.posterior.data_vars:
                    init_values_chain = dict()
                    init_values_chain[var_name] = vi_samples.posterior[var_name].to_numpy()
                init_values.append(init_values_chain)

            trace = pm.sample(
                draws=samples,
                tune=tunings,
                chains=n_chains,
                cores=1,
                target_accept=target_accept,
                initvals=init_values,
                return_inferencedata=True,
                progressbar=False
            )
        n_diverging = sum(sum(trace.sample_stats.diverging.values))
        if not target_accept_low:
            non_div_run = True
        elif n_diverging > 0:
            target_accept = .99
            tunings += 500
            target_accept_low = False
        else:
            non_div_run = True

    beta_detection_mu = trace.posterior["beta_detection"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["beta_detection"].std(("chain", "draw")).values

    ordered_values = [val for pair in zip(beta_detection_mu, beta_detection_sd) for val in pair]

    beta_bait_mu = trace.posterior["beta_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["beta_bait"].std(("chain", "draw")).item()
    credible_interval = az.hdi(trace, var_names=["beta_bait"], hdi_prob=0.95)
    low_ci, hugh_ci = credible_interval["beta_bait"].values

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    detection_parameters = ordered_values + [beta_bait_mu, beta_bait_sd, low_ci, hugh_ci, n_diverging]
    end = datetime.now()
    print(f"A job took {end - start} time at target_accept {target_accept} of job{j}")
    # Return a tab-separated string for writing later
    return "\t".join(map(str, interaction_row.values)) + "\t" + "\t".join(map(str, detection_parameters))


@ray.remote(num_cpus=1)
def process_prey_abundance(interaction_row, in_obs_c, in_tested_c, cl_categories, untargeted_df, j, samples=1500,
                           tunings=1000, smallstep=False):
    import pymc as pm
    import logging
    logger = logging.getLogger("pymc")
    logger.setLevel(logging.ERROR)
    from pymc.variational.callbacks import CheckParametersConvergence


    mask = np.zeros(cl_categories.shape, dtype=bool)
    for i, cl_tested in enumerate(in_tested_c):
        if interaction_row[cl_tested] != 0:
            mask[i] = True

    fitted_cl_categories = cl_categories[mask]
    n_categories = int(sum(mask))
    tested_c = in_tested_c[mask]
    obs_c = in_obs_c[mask]

    start = datetime.now()
    if smallstep:
        target_accept = 0.99
    else:
        target_accept = 0.9
    n_tests = interaction_row[tested_c].sum()
    bait_matrix = np.zeros((n_tests, n_categories + 2))
    bait_matrix[:, 1] = 1
    i = 0
    cl_i = 2
    cl_idx = np.zeros(bait_matrix.shape[0], dtype=int)

    for cl_obs, cl_test in zip(obs_c, tested_c):
        cl_idx[i:(i + interaction_row[cl_test])] = cl_i - 2
        bait_matrix[i:(i + interaction_row[cl_obs]), 0] = 1
        bait_matrix[i:(i + interaction_row[cl_test]), cl_i] = 1
        i += interaction_row[cl_test]
        cl_i += 1

    untargeted_df = untargeted_df[untargeted_df["cell_line"].isin(fitted_cl_categories)].copy()
    category_to_index = {cat: idx for idx, cat in enumerate(fitted_cl_categories)}
    untargeted_df['cell_line_idx'] = untargeted_df['cell_line'].map(category_to_index)
    observed_x = untargeted_df.iloc[:, 0].values
    group_idx = untargeted_df['cell_line_idx'].values

    target_accept_low = True
    non_div_run = False
    n_chains = 4

    while not non_div_run:
        with pm.Model() as multi_env_model:
            mu_prior = 0
            sd_prior = 1  # as we harmonised

            mu = pm.Normal('mu', mu=mu_prior, sigma=sd_prior, shape=n_categories)
            sigma = pm.HalfNormal('sigma', sigma=sd_prior, shape=n_categories)  # or another prior on std dev

            x_untargeted = pm.Normal(
                'x_untargeted',
                mu=mu[group_idx],
                sigma=sigma[group_idx],
                observed=observed_x)

            x_to_samples = x_untargeted[cl_idx]
            b_bait = pm.SkewNormal('b_bait', mu=0, sigma=10)

            p_y = pm.Deterministic("mu_y", pm.math.sigmoid(b_bait * (2 ** x_to_samples)))  # we estimate mean on log
            y_lh = pm.Bernoulli('y_lh', p=p_y, observed=bait_matrix[:, 0])

            vi_approx = pm.fit(
                method="advi",
                callbacks=[CheckParametersConvergence(diff="absolute")])
            init_values = []
            for _ in range(n_chains):
                vi_samples = vi_approx.sample(1)
                for var_name in vi_samples.posterior.data_vars:
                    init_values_chain = dict()
                    init_values_chain[var_name] = vi_samples.posterior[var_name].to_numpy()
                init_values.append(init_values_chain)

            trace = pm.sample(
                draws=samples,
                tune=tunings,
                chains=n_chains,
                cores=1,
                target_accept=target_accept,
                initvals=init_values,
                return_inferencedata=True,
                progressbar=False
            )
        n_diverging = sum(sum(trace.sample_stats.diverging.values))
        if not target_accept_low:
            non_div_run = True
        elif n_diverging > 0:
            target_accept = .99
            tunings += 500
            target_accept_low = False
        else:
            non_div_run = True

    beta_detection_mu = trace.posterior["mu"].mean(("chain", "draw")).values
    beta_detection_sd = trace.posterior["sigma"].std(("chain", "draw")).values

    ordered_values = []
    current = 0
    for included in mask:
        if included:
            ordered_values.append(beta_detection_mu[current])
            ordered_values.append(beta_detection_sd[current])
            current += 1
        else:
            ordered_values += [np.nan] * 2

    beta_bait_mu = trace.posterior["b_bait"].mean(("chain", "draw")).item()
    beta_bait_sd = trace.posterior["b_bait"].std(("chain", "draw")).item()

    credible_interval = az.hdi(trace, var_names=["b_bait"], hdi_prob=0.95)
    low_ci, hugh_ci = credible_interval["b_bait"].values

    n_diverging = sum(sum(trace.sample_stats.diverging.values))

    detection_parameters = ordered_values + [beta_bait_mu, beta_bait_sd, low_ci, hugh_ci, n_diverging]
    end = datetime.now()
    print(f"A job took {end - start} time at target_accept {target_accept} of job{j}")
    # Return a tab-separated string for writing later
    return "\t".join(map(str, interaction_row.values)) + "\t" + "\t".join(map(str, detection_parameters))


def main():
    parser = argparse.ArgumentParser(description="Process bait interaction data with multiprocessing.")
    parser.add_argument("--prey_tested", required=True, help="Path to prey combination TSV file")
    parser.add_argument("--abundance_cell_lines", required=True, help="Path to baseline pod TSV file")
    parser.add_argument("--bait_output", required=True, help="Path to output file for bait parameters")
    parser.add_argument("--abundance", type=int, help="Set to 1 for abundance, 0 for probability of detection")
    parser.add_argument("--samples", type=int, help="Number of samples from posterior")
    parser.add_argument("--burin_samples", type=int, help="Number of burnin samples for posterior")
    parser.add_argument("--workers", type=int, help="Number of worker processes for ray")
    parser.add_argument("--batch_size", type=int, help="Number of tests performed before writing to file")
    parser.add_argument("--stepsize", type=argparse.BooleanOptionalAction, help="sets minimal stepsize")
    args = parser.parse_args()
    smallstep = args.stepsize

    prey_interaction_df = pd.read_csv(args.prey_tested, sep="\t")
    numeric_cols = prey_interaction_df.columns[1:]
    prey_interaction_df[numeric_cols] = prey_interaction_df[numeric_cols].astype(int)
    cell_line_abundance = pd.read_csv(args.abundance_cell_lines, sep="\t")
    cl_categories = np.array(["CVCL_0030", "CVCL_0063", "CVCL_0291"])  # TODO set as argument

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
                prey_abundance_df = cell_line_abundance[[prey, "cell_line"]].dropna().copy()
                futures.append(process_prey_abundance.options(
                    runtime_env={"env_vars": ray_task_env}
                ).remote(
                    interaction_row=row,
                    in_obs_c=observed_cols,
                    in_tested_c=tested_cols,
                    cl_categories=cl_categories,
                    untargeted_df=prey_abundance_df,
                    j=j,
                    samples=args.samples,
                    tunings=args.burin_samples,
                    smallstep=smallstep
                ))

            elif args.abundance == 0:
                detection_matrix = get_one_hot_untargeted(
                    cl_categories,
                    cell_line_abundance[[prey, "cell_line"]].fillna(0))  # remove non-observations

                futures.append(process_prey_pod.options(
                    runtime_env={"env_vars": ray_task_env}
                ).remote(
                    interaction_row=row,
                    obs_c=observed_cols,
                    tested_c=tested_cols,
                    cl_categories=cl_categories,
                    detection_matrix=detection_matrix,
                    j=j,
                    samples=args.samples,
                    tunings=args.burin_samples
                ))

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
