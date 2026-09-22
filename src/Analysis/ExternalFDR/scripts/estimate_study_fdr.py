import os
os.environ["OMP_NUM_THREADS"] = "1"   # one thread per job, the pool is the parallelism

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize
import multiprocessing as mp

def fit_mixture(bg_df, max_iter, tol=1e-8):
    # Run across unique observed/tested
    cell = bg_df.groupby(["n_observed", "n_tested"]).size().reset_index(name="w")
    n, t, w = cell.n_observed.values, cell.n_tested.values, cell.w.values.astype(float) # w is weight as it the number of pairs with same ratio

    pi, a, b, phi = 0.01, 0.5, 2.5, 4e-4 # EM starting values, pi: true interaction rate, a: alpha value, b: beta value, phi:FPR
    ll_old = -np.inf

    for _ in range(max_iter):
        A = betabinom.pmf(n, t, a, b)
        B = binom.pmf(n, t, phi)
        den = pi * A + (1 - pi) * B
        g = pi * A / den                                    # E-step

        wg, wng = w * g, w * (1 - g)                        # M-step
        pi  = wg.sum() / w.sum()
        phi = (wng * n).sum() / (wng * t).sum()

        nll = lambda x: -(wg * betabinom.logpmf(n, t, *np.exp(x))).sum()
        a, b = np.exp(minimize(nll, np.log([a, b]), method="Nelder-Mead").x)

        ll = (w * np.log(pi * betabinom.pmf(n, t, a, b)
                         + (1 - pi) * binom.pmf(n, t, phi))).sum()  # den is pre M-step
        if ll - ll_old < tol * abs(ll):                             # relative, loglik ~1e8
            ll_old = ll
            break
        ll_old = ll

    return dict(pi=pi, a=a, b=b, phi=phi,
                p_mean=a / (a + b),
                p_var=a * b / ((a + b) ** 2 * (a + b + 1)),
                loglik=ll_old)
    

def score_interactions_single_study(external_study_df, g_pi, g_a, g_b, g_phi, max_iter, tol=1e-8):
    n, t = external_study_df.n_observed.values, external_study_df.n_tested.values

    A = betabinom.pmf(n, t, g_a, g_b)      # external evidence | true
    B = binom.pmf(n, t, g_phi)             # external evidence | false positive

    theta = g_pi                           # start at the population prevalence
    for _ in range(max_iter):
        g = theta * A / (theta * A + (1 - theta) * B)
        new = g.mean()
        if abs(new - theta) < tol:
            theta = new
            break
        theta = new

    L = theta * A + (1 - theta) * B
    se = 1.0 / np.sqrt(np.sum(((A - B) / L) ** 2))

    return dict(study_theta=theta, fdr=1 - theta, se=se,
                n_external_reported=len(n))
    

def _run_study(study):
    study_name = os.path.basename(study).removesuffix(".csv")
    s_df = pd.read_csv(study, sep = "\t")
    if not directed:
        s_df["pair_id"] = s_df.apply(lambda row: 
            id_dict.get(tuple(sorted([row[prot_a], row[prot_b]])),
                        -1), axis=1)
        
        s_df = s_df.groupby("pair_id", as_index=False, sort=False).agg({
                "n_tested": "sum",
                "n_observed": "sum"
            })
    else:
        s_df["pair_id"] = s_df.apply(lambda row: 
            id_dict.get((row[prot_a], row[prot_b]), -1), axis=1)
    
    s_df = s_df[s_df.pair_id != -1]
    s_df.set_index("pair_id", inplace=True)
    detected_pairs = s_df[s_df.n_observed > 0].index 
    
    background_df = full_search_df.copy()
    background_df.loc[s_df.index, "n_tested"] -= s_df.n_tested
    background_df.loc[detected_pairs, "n_observed"] -= s_df.loc[detected_pairs, "n_observed"]
    background_df = background_df[background_df.n_tested != 0]
    
    assert (background_df.n_tested < 1).sum() == 0, "There are negative tested"
    
    external_detected = background_df[background_df.index.isin(detected_pairs)]
    
    global_estimates = fit_mixture(background_df, max_iter)
    study_estimates = score_interactions_single_study(
        external_detected,
        global_estimates["pi"],
        global_estimates["a"],
        global_estimates["b"],
        global_estimates["phi"],
        max_iter)
    
    row_dict = {"study": study_name}
    row_dict.update(study_estimates)
    row_dict.update(global_estimates)
    return row_dict


def main():
    global full_search_df, id_dict, directed, prot_a, prot_b, max_iter
    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    full_search_df = pd.read_parquet(
        snakemake.input.pod_file,
        columns=[prot_a, prot_b, "pair_id", "n_tested", "n_observed"])
    consituent_studies = snakemake.input.consituent_studies

    max_iter = snakemake.params.n_em_iterations
    n_threads = snakemake.threads

    if snakemake.wildcards.network_type=="undirectional":
        directed=False
    elif snakemake.wildcards.network_type=="directional":
        directed=True
    else:
        raise IOError(f"Unkown network wildcard: {snakemake.wildcards.network_type}")

    id_dict = {(row[prot_a], row[prot_b]): row["pair_id"] for _, row in full_search_df.iterrows()}
    # only the counts are read per study, drop the rest before the workers copy it
    full_search_df = full_search_df.set_index("pair_id")[["n_tested", "n_observed"]]

    with mp.Pool(n_threads) as pool:
        study_rows = pool.map(_run_study, consituent_studies, chunksize=1)

    study_fdrs = pd.DataFrame(study_rows)
    study_fdrs.to_csv(snakemake.output.study_fdr, sep="\t", index=False)
    
if __name__== "__main__":
    main()