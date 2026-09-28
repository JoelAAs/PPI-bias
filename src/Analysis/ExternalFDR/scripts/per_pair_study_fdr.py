import os

os.environ["OMP_NUM_THREADS"] = "1"  # one thread per job, the pool is the parallelism

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize
import multiprocessing as mp

# ----------------------------------------------------------------------------
# Phase 0: global mixture fit
# ----------------------------------------------------------------------------


def combinations_count(n_observed, n_tested):
    """
    Return unique key per ntested/observed combinations and the number of studies with that combinations (for weighting)
    """
    keys, weights = np.unique(n_observed * max_tested + n_tested, return_counts=True)
    return pd.Series(weights, index=keys)


def fit_mixture(pair_testing_combinations, max_iter, tol=1e-8):
    n, t = pair_testing_combinations.index.values // max_tested, pair_testing_combinations.index.values % max_tested
    w = pair_testing_combinations.values.astype(float)

    pi, a, b, phi = 0.01, 0.5, 2.5, 4e-4
    ll_old = -np.inf
    converged = False

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
            converged = True
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
        # true pairs must be detected more often than false ones, else the components swapped
        swapped=a / (a + b) <= phi,
        converged=converged,  # False when max_iter was hit before tol
    )



def score_reported(n, t, g_pi, g_a, g_b, g_phi, max_iter, tol=1e-8):
    """
    Given the background estimation of alpha beta, interaction-rate and false negative rate
    Given whats reported in this study, what the precision of the study?
    """
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


def score_negatives(n, t, g_pi, g_a, g_b, g_phi):
    """
    From the rest of the studies, what is the expected number of false negatives.
    """
    A = betabinom.pmf(n, t, g_a, g_b) # P(external evidence | pair is a true interaction)
    B = binom.pmf(n, t, g_phi)        # P(external evidence | pair is not an interaction)
    gamma = g_pi * A / (g_pi * A + (1 - g_pi) * B) # P(true interaction | external evidence)
    return dict(for_g_sum=gamma.sum(), n_negative_informative=len(n))


def lookup_pair_ids(bait, prey):
    """ protein ids -> POD pair_id via the sorted pair-key table, -1 when not in POD """
    ca, cb = id_vocab.get_indexer(bait), id_vocab.get_indexer(prey)
    if not directed:
        # vocab is sorted, so ordering codes is ordering ids
        ca, cb = np.minimum(ca, cb), np.maximum(ca, cb)
    key = ca.astype(np.int64) * len(id_vocab) + cb
    pos = np.searchsorted(pod_keys, key).clip(max=len(pod_keys) - 1)
    hit = (ca >= 0) & (cb >= 0) & (pod_keys[pos] == key)
    return np.where(hit, pod_key_pair_ids[pos], -1)


def load_study(study_path):
    s_df = pd.read_csv(
        study_path, sep="\t", usecols=[prot_a, prot_b, "n_tested", "n_observed"]
    )
    s_df["pair_id"] = lookup_pair_ids(s_df[prot_a].values, s_df[prot_b].values)
    s_df = s_df[s_df.pair_id != -1]
    s_df = s_df.groupby("pair_id", sort=False)[["n_tested", "n_observed"]].sum()
    return (
        s_df.index.values.astype(np.int32),
        s_df.n_observed.values.astype(np.int32),
        s_df.n_tested.values.astype(np.int32),
    )


def _run_study(study):
    study_name = os.path.basename(study).removesuffix(".csv")
    pair_ids, s_obs, s_tst = load_study(study)
    detected = s_obs > 0

    # Leave-one-out background, only the pairs this study touched change
    old_o, old_t = pod_observed[pair_ids], pod_tested[pair_ids]
    new_o, new_t = old_o - s_obs, old_t - s_tst
    assert (new_t < 0).sum() == 0, "There are negative tested"
    kept = new_t != 0

    background = (
        bg_cells.sub(combinations_count(old_o, old_t), fill_value=0)
        .add(combinations_count(new_o[kept], new_t[kept]), fill_value=0)
    )
    assert (background < 0).sum() == 0, "Negative cell count in background"
    background = background[background > 0]

    global_estimates = fit_mixture(background, max_iter)
    gp = (
        global_estimates["pi"],
        global_estimates["a"],
        global_estimates["b"],
        global_estimates["phi"],
    )

    row = dict(
        study=study_name,
        k_s=float(s_obs.sum()),  # positive calls made by this study
        T_s=float(s_tst.sum()),  # tests performed by this study
        n_reported=int(detected.sum()),
        n_negative=int((~detected).sum()),
    )
    # Interactions reported in backgound, for FDR
    rep = detected & kept
    row.update(score_reported(new_o[rep], new_t[rep], *gp, max_iter))
    # Tested but not reported pairs for FOR
    neg = ~detected & kept
    row.update(score_negatives(new_o[neg], new_t[neg], *gp))
    row.update(global_estimates)

    # keep the study's counts for phase 2 so it is not loaded twice
    return row, (pair_ids, s_obs, s_tst)


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


def accumulate_log_odds(study_counts, weights, n_pairs, prior_logodds):
    total = np.full(n_pairs, prior_logodds, dtype=np.float64)
    for study, (pair_ids, obs, tst) in study_counts.items():
        w_det, w_non = weights.loc[study, "w_detected"], weights.loc[study, "w_not_detected"]
        if np.isnan(w_det):
            continue
        # each test is an independent opportunity; detections carry w_det, the rest w_non
        # pair_ids are unique within a study, so plain fancy-index add is safe
        total[pair_ids] += obs * w_det + (tst - obs) * w_non
    return total


def main():
    global directed, prot_a, prot_b, max_iter, max_tested
    global id_vocab, pod_keys, pod_key_pair_ids, pod_observed, pod_tested
    global bg_cells
    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    pod_df = pd.read_parquet(
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

    # (bait, prey) -> pair_id as a sorted int64 key table instead of a dict of tuples:
    # numpy arrays stay shared across forked workers, a 1e8 entry dict gets copied
    n_rows = len(pod_df)
    codes, id_vocab = pd.factorize(
        np.concatenate([pod_df[prot_a].values, pod_df[prot_b].values]), sort=True
    )
    id_vocab = pd.Index(id_vocab)
    keys = codes[:n_rows].astype(np.int64) * len(id_vocab) + codes[n_rows:]
    order = np.argsort(keys)
    pod_keys = keys[order]
    pod_key_pair_ids = pod_df.pair_id.values[order]
    del codes, keys, order

    # counts as dense arrays indexed by pair_id
    pod_pair_ids = pod_df.pair_id.values
    n_pairs = int(pod_pair_ids.max()) + 1
    assert np.bincount(pod_pair_ids, minlength=n_pairs).max() <= 1, "Duplicate pair_id in POD"
    pod_tested = np.zeros(n_pairs, dtype=np.int64)
    pod_observed = np.zeros(n_pairs, dtype=np.int64)
    pod_tested[pod_pair_ids] = pod_df.n_tested.values
    pod_observed[pod_pair_ids] = pod_df.n_observed.values
    del pod_df

    max_tested = int(pod_tested.max()) + 1
    t_all, o_all = pod_tested[pod_pair_ids], pod_observed[pod_pair_ids]

    # global prior from the full fit
    global_fit = fit_mixture(combinations_count(o_all, t_all), max_iter)
    bg_cells = combinations_count(o_all[t_all != 0], t_all[t_all != 0])
    del t_all, o_all

    # ---- phase 1: per-study FDR and FOR --------------------------------------
    study_rows, study_counts = [], {}
    with mp.Pool(n_threads) as pool:
        for row, counts in pool.imap(_run_study, consituent_studies, chunksize=1):
            study_rows.append(row)
            study_counts[row["study"]] = counts

    study_df = pd.DataFrame(study_rows)
    study_df = solve_rates(study_df, global_fit["pi"], k0=k0)
    study_df.to_csv(snakemake.output.study_fdr, sep="\t", index=False)

    # ---- phase 2: per-pair posterior ----------------------------------------
    weights = study_df.set_index("study")[["w_detected", "w_not_detected"]]

    pi0 = global_fit["pi"]
    prior_lo = np.log(pi0 / (1 - pi0))

    log_odds = accumulate_log_odds(study_counts, weights, n_pairs, prior_lo)

    pairs = pd.DataFrame(
        {
            "pair_id": pod_pair_ids,
            "n_tested": pod_tested[pod_pair_ids],
            "n_observed": pod_observed[pod_pair_ids],
            "log_odds": log_odds[pod_pair_ids],
        }
    )
    pairs["posterior"] = 1.0 / (1.0 + np.exp(-pairs.log_odds))
    pairs.to_parquet(snakemake.output.pair_posterior, index=False)


if __name__ == "__main__":
    main()
