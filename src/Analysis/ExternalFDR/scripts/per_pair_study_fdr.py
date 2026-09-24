import os

os.environ["OMP_NUM_THREADS"] = "1"  # one thread per job, the pool is the parallelism

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize
import multiprocessing as mp

# ----------------------------------------------------------------------------
# Phase 0: global mixture fit  (unchanged)
# ----------------------------------------------------------------------------


def fit_mixture(bg_df, max_iter, tol=1e-8):
    cell = bg_df.groupby(["n_observed", "n_tested"]).size().reset_index(name="w")
    n, t, w = cell.n_observed.values, cell.n_tested.values, cell.w.values.astype(float)

    pi, a, b, phi = 0.01, 0.5, 2.5, 4e-4
    ll_old = -np.inf

    for _ in range(max_iter):
        A = betabinom.pmf(n, t, a, b)
        B = binom.pmf(n, t, phi)
        den = pi * A + (1 - pi) * B
        g = pi * A / den  # E-step

        wg, wng = w * g, w * (1 - g)  # M-step
        pi = wg.sum() / w.sum()
        phi = (wng * n).sum() / (wng * t).sum()

        nll = lambda x: -(wg * betabinom.logpmf(n, t, *np.exp(x))).sum()
        a, b = np.exp(minimize(nll, np.log([a, b]), method="Nelder-Mead").x)

        ll = (
            w * np.log(pi * betabinom.pmf(n, t, a, b) + (1 - pi) * binom.pmf(n, t, phi))
        ).sum()
        if ll - ll_old < tol * abs(ll):
            ll_old = ll
            break
        ll_old = ll

    return dict(
        pi=pi,
        a=a,
        b=b,
        phi=phi,
        p_mean=a / (a + b),
        p_var=a * b / ((a + b) ** 2 * (a + b + 1)),
        loglik=ll_old,
    )



def score_reported(external_df, g_pi, g_a, g_b, g_phi, max_iter, tol=1e-8):
    """
    Given the background estimation of alpha beta, interaction-rate and false negative rate
    Given whats reported in this study, what the precision of the study?
    """
    n, t = external_df.n_observed.values, external_df.n_tested.values
    if len(n) == 0:
        return dict(theta_g_sum=np.nan, n_reported_informative=0, theta_se=np.inf)

    A = betabinom.pmf(n, t, g_a, g_b)  # external evidence | true
    B = binom.pmf(n, t, g_phi)  # external evidence | false positive

    theta = g_pi
    for _ in range(max_iter):
        g = theta * A / (theta * A + (1 - theta) * B) # True prob / reported true prob
        new = g.mean()
        if abs(new - theta) < tol:
            theta = new
            break
        theta = new

    L = theta * A + (1 - theta) * B
    with np.errstate(divide="ignore"):
        se = 1.0 / np.sqrt(np.sum(((A - B) / L) ** 2))

    # return the sufficient statistic, not the point estimate: main() shrinks it
    return dict(theta_g_sum=g.sum(), n_reported_informative=len(n), theta_se=se)


def score_negatives(external_df, g_pi, g_a, g_b, g_phi):
    """
    From the rest of the studies, what is the expected number of false negatives. 
    """
    n, t = external_df.n_observed.values, external_df.n_tested.values
    
    A = betabinom.pmf(n, t, g_a, g_b) # P(external evidence | pair is a true interaction)
    B = binom.pmf(n, t, g_phi)        # P(external evidence | pair is not an interaction)
    gamma = g_pi * A / (g_pi * A + (1 - g_pi) * B) # P(true interaction | external evidence)
    return dict(for_g_sum=gamma.sum(), n_negative_informative=len(n))


def load_study(study_path):
    s_df = pd.read_csv(study_path, sep="\t")
    if not directed:
        s_df["pair_id"] = s_df.apply(
            lambda row: id_dict.get(tuple(sorted([row[prot_a], row[prot_b]])), -1),
            axis=1,
        )
        s_df = s_df.groupby("pair_id", as_index=False, sort=False).agg(
            {"n_tested": "sum", "n_observed": "sum"}
        )
    else:
        s_df["pair_id"] = s_df.apply(
            lambda row: id_dict.get((row[prot_a], row[prot_b]), -1), axis=1
        )

    s_df = s_df[s_df.pair_id != -1]
    return s_df.set_index("pair_id")


def _run_study(study):
    study_name = os.path.basename(study).removesuffix(".csv")
    s_df = load_study(study)

    detected_pairs = s_df.index[s_df.n_observed > 0]
    k_s = float(s_df.n_observed.sum())  # positive calls made by this study
    T_s = float(s_df.n_tested.sum())  # tests performed by this study

    # Leave-one-out background
    background_df = full_search_df.copy()
    background_df.loc[s_df.index, "n_tested"] -= s_df.n_tested
    background_df.loc[detected_pairs, "n_observed"] -= s_df.loc[
        detected_pairs, "n_observed"
    ]
    background_df = background_df[background_df.n_tested != 0]

    assert (background_df.n_tested < 1).sum() == 0, "There are negative tested"

    global_estimates = fit_mixture(background_df, max_iter)
    gp = (
        global_estimates["pi"],
        global_estimates["a"],
        global_estimates["b"],
        global_estimates["phi"],
    )

    # Interactions reported in backgound, for FDR
    reported_ext = background_df.loc[detected_pairs]

    # Tested but not reported pairs for FOR
    negative_pairs = s_df.index[s_df.n_observed == 0]
    negative_ext = full_search_df.loc[negative_pairs].copy()

    negative_ext["n_tested"] -= s_df.loc[negative_ext.index, "n_tested"]
    negative_ext = negative_ext[negative_ext.n_tested != 0]

    row = dict(
        study=study_name,
        k_s=k_s,
        T_s=T_s,
        n_reported=len(detected_pairs),
        n_negative=len(negative_pairs),
    )
    row.update(score_reported(reported_ext, *gp, max_iter))
    row.update(score_negatives(negative_ext, *gp))
    row.update(global_estimates)
    return row


# ----------------------------------------------------------------------------
# Phase 1 post-processing: shrinkage, then FDR/FOR -> p_s, phi_s
# ----------------------------------------------------------------------------


def shrink(g_sum, n, prior_mean, k0):
    """  we shift towards mean in order to not get infinete OR (low overlap studies) """
    return (g_sum + k0 * prior_mean) / (n + k0)


def solve_rates(df, global_pi, k0=5.0, eps=1e-9):
    
    # Mean precision weighted by number of overlap from estimation
    m_theta = np.nansum(df.theta_g_sum) / np.nansum(df.n_reported_informative)
    # Mean FOR weighted by number of overlap from estimation
    m_for = np.nansum(df.for_g_sum) / np.nansum(df.n_negative_informative)

    theta = shrink(df.theta_g_sum.fillna(0), df.n_reported_informative, m_theta, k0)
    FOR = shrink(df.for_g_sum.fillna(0), df.n_negative_informative, m_for, k0)

    r = df.k_s / df.T_s # call rate
    FDR = 1 - theta

    pi_s = r * theta + FOR * (1 - r)
    pi_s = pi_s.clip(eps, 1 - eps)
    p_s = (r * theta / pi_s).clip(eps, 1 - eps) # sensitivity
    phi_s = (r * FDR / (1 - pi_s)).clip(eps, 1 - eps)

    # a study cannot be less sensitive than it is spurious; flag rather than hide
    bad = (p_s >= 1 - 1e-6) | p_s.isna()
    p_s = np.where(bad, np.nan, p_s)
    phi_s = np.where(bad, np.nan, phi_s)

    out = df.copy()
    out["theta_s"], out["fdr_s"], out["for_s"] = theta, FDR, FOR
    out["r_s"], out["pi_s"], out["p_s"], out["phi_s"] = r, pi_s, p_s, phi_s
    out["w_detected"] = np.log(p_s / phi_s)
    out["w_not_detected"] = np.log((1 - p_s) / (1 - phi_s))
    out["degenerate"] = bad
    return out


def _pair_contrib(args):
    study, w_det, w_non = args
    s_df = load_study(study)
    obs = s_df.n_observed.values
    tst = s_df.n_tested.values
    # each test is an independent opportunity; detections carry w_det, the rest w_non
    contrib = obs * w_det + (tst - obs) * w_non
    return s_df.index.values, contrib


def accumulate_log_odds(study_paths, weights, n_pairs, prior_logodds, n_threads):
    total = np.full(n_pairs, prior_logodds, dtype=np.float64)
    jobs = [
        (s, weights.loc[s, "w_detected"], weights.loc[s, "w_not_detected"])
        for s in study_paths
        if not np.isnan(weights.loc[s, "w_detected"])
    ]

    with mp.Pool(n_threads) as pool:
        for ids, contrib in pool.imap_unordered(_pair_contrib, jobs, chunksize=1):
            np.add.at(total, ids, contrib)
    return total


def main():
    global full_search_df, id_dict, directed, prot_a, prot_b, max_iter
    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    full_search_df = pd.read_parquet(
        snakemake.input.pod_file,
        columns=[prot_a, prot_b, "pair_id", "n_tested", "n_observed"],
    )
    consituent_studies = snakemake.input.consituent_studies

    max_iter = snakemake.params.n_em_iterations
    n_threads = snakemake.threads
    k0 = getattr(snakemake.params, "shrinkage_pseudocount", 5.0)

    if snakemake.wildcards.network_type == "undirectional":
        directed = False
    elif snakemake.wildcards.network_type == "directional":
        directed = True
    else:
        raise IOError(f"Unkown network wildcard: {snakemake.wildcards.network_type}")

    id_dict = {
        (row[prot_a], row[prot_b]): row["pair_id"]
        for _, row in full_search_df.iterrows()
    }
    full_search_df = full_search_df.set_index("pair_id")[["n_tested", "n_observed"]]

    # ---- phase 1: per-study FDR and FOR --------------------------------------
    with mp.Pool(n_threads) as pool:
        study_rows = pool.map(_run_study, consituent_studies, chunksize=1)

    study_df = pd.DataFrame(study_rows)

    # global prior from the full fit
    global_fit = fit_mixture(full_search_df, max_iter)
    study_df = solve_rates(study_df, global_fit["pi"], k0=k0)
    study_df.to_csv(snakemake.output.study_fdr, sep="\t", index=False)

    # ---- phase 2: per-pair posterior ----------------------------------------
    weights = study_df.set_index("study")[["w_detected", "w_not_detected"]]
    weights.index = [
        os.path.join(os.path.dirname(consituent_studies[0]), f"{s}.csv")
        for s in weights.index
    ]

    pi0 = global_fit["pi"]
    prior_lo = np.log(pi0 / (1 - pi0))
    n_pairs = int(full_search_df.index.max()) + 1

    log_odds = accumulate_log_odds(
        consituent_studies, weights, n_pairs, prior_lo, n_threads
    )

    pairs = pd.DataFrame(
        {
            "pair_id": full_search_df.index.values,
            "n_tested": full_search_df.n_tested.values,
            "n_observed": full_search_df.n_observed.values,
            "log_odds": log_odds[full_search_df.index.values],
        }
    )
    pairs["posterior"] = 1.0 / (1.0 + np.exp(-pairs.log_odds))
    pairs.to_parquet(snakemake.output.pair_posterior, index=False)


if __name__ == "__main__":
    main()
