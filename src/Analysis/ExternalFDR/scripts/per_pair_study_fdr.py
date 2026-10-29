import os

os.environ["OMP_NUM_THREADS"] = "1"  # one thread per job, the pool is the parallelism

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logit
import multiprocessing as mp


def unique_test_comb(n_observed, n_tested):
    """(n_observed, n_tested) -> unique combinations and how many pairs hold each."""
    keys, w = np.unique(n_observed * max_tested + n_tested, return_counts=True)
    return keys // max_tested, keys % max_tested, w.astype(float)


def fit_mixture(n, t, w, max_iter, tol=1e-8):
    """Two components: true pairs ~ BetaBinom(t, a, b), false positives ~ Binom(t, phi)."""
    pi, a, b, phi = 0.01, 0.5, 2.5, 4e-4
    ll = ll_old = -np.inf
    converged = False

    for _ in range(max_iter):
        A, B = betabinom.pmf(n, t, a, b), binom.pmf(n, t, phi)
        den = pi * A + (1 - pi) * B
        g = pi * A / den  # E-step

        wg, wng = w * g, w * (1 - g)  # M-step
        pi = wg.sum() / w.sum()
        phi = (wng * n).sum() / (wng * t).sum()
        nll = lambda x: -(wg * betabinom.logpmf(n, t, *np.exp(x))).sum()
        a, b = np.exp(minimize(nll, np.log([a, b]), method="Nelder-Mead").x)

        ll = (w * np.log(den)).sum()
        if ll - ll_old < tol * abs(ll):
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
        loglik=ll,
        swapped=a / (a + b) <= phi,  # else the components changed places
        converged=converged,
    )


T_BINS = np.array([1, 2, 3, 4, 5])  # lower edges of the external-test bins


def confirm_counts(n, t):
    """per bin: [informative pairs, pairs with >= 1 external hit]"""
    b = np.digitize(t, T_BINS) - 1
    return np.stack(
        [
            np.bincount(b, minlength=len(T_BINS)),
            np.bincount(b, n > 0, minlength=len(T_BINS)),
        ]
    )


def score_reported(n, t, c, phi):
    """theta_s by ML: an informative pair is confirmed (>= 1 external hit) with
    probability c(t) if true and 1 - (1 - phi)^t if a false positive."""
    if len(n) == 0:
        return dict(theta_g_sum=np.nan, n_informative=0, mean_t_ext=np.nan)
    
    conf = n > 0
    q = 1 - (1 - phi) ** t # chance of at least one FP
    ct = np.clip(c[np.digitize(t, T_BINS) - 1], q + 1e-12, 1 - 1e-12)
    likelihood = lambda theata: np.where(conf, theata * ct + (1 - th) * q, 1 - theata * ct - (1 - th) * q)
    theta = minimize_scalar(
        lambda th: -np.log(likelihood(th)).sum(), bounds=(0, 1), method="bounded"
    ).x
    g = theta * np.where(conf, ct, 1 - ct) / lik(theta)  # P(true | external evidence)
    # the sufficient statistic, not the point estimate: main() shrinks it
    return dict(
        theta_g_sum=g.sum(),
        n_informative=len(n),
        mean_t_ext=t.mean(),
        f_confirmed=conf.mean(),
        f_expected_true=ct.mean(),
    )


def lookup_pair_ids(bait, prey):
    """protein ids -> POD pair_id via the sorted pair-key table, -1 when not in POD."""
    ca, cb = id_vocab.get_indexer(bait), id_vocab.get_indexer(prey)
    if not directed:
        ca, cb = np.minimum(ca, cb), np.maximum(
            ca, cb
        )  # vocab sorted: codes order as ids
    key = ca.astype(np.int64) * len(id_vocab) + cb
    pos = np.searchsorted(pod_keys, key).clip(max=len(pod_keys) - 1)
    hit = (ca >= 0) & (cb >= 0) & (pod_keys[pos] == key)
    return np.where(hit, pod_key_pair_ids[pos], -1)


def load_study(path):
    df = pd.read_csv(path, sep="\t", usecols=[prot_a, prot_b, "n_tested", "n_observed"])
    df["pair_id"] = lookup_pair_ids(df[prot_a].values, df[prot_b].values)
    df = (
        df[df.pair_id != -1]
        .groupby("pair_id", sort=False)[["n_tested", "n_observed"]]
        .sum()
    )
    return (
        os.path.basename(path).removesuffix(".csv"),
        df.index.values.astype(np.int32),
        df.n_observed.values.astype(np.int32),
        df.n_tested.values.astype(np.int32),
    )


def main():
    global directed, prot_a, prot_b, max_tested
    global id_vocab, pod_keys, pod_key_pair_ids

    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    max_iter = snakemake.params.n_em_iterations
    k0 = getattr(snakemake.params, "shrinkage_pseudocount", 5.0)
    n_ref = getattr(snakemake.params, "n_reference_studies", 10)
    directed = {"undirectional": False, "directional": True}[
        snakemake.wildcards.network_type
    ]

    pod_df = pd.read_parquet(
        snakemake.input.pod_file,
        columns=[prot_a, prot_b, "pair_id", "n_tested", "n_observed"],
    )

    # (bait, prey) -> pair_id as a sorted int64 key table: numpy arrays stay shared
    # across forked workers, a 1e8 entry dict gets copied
    n_rows = len(pod_df)
    codes, id_vocab = pd.factorize(
        np.concatenate([pod_df[prot_a].values, pod_df[prot_b].values]), sort=True
    )
    id_vocab = pd.Index(id_vocab)
    keys = codes[:n_rows].astype(np.int64) * len(id_vocab) + codes[n_rows:]
    order = np.argsort(keys)
    pod_keys, pod_key_pair_ids = keys[order], pod_df.pair_id.values[order]

    # counts as dense arrays indexed by pair_id
    pair_ids = pod_df.pair_id.values
    n_pairs = int(pair_ids.max()) + 1
    assert (
        np.bincount(pair_ids, minlength=n_pairs).max() <= 1
    ), "Duplicate pair_id in POD"
    pod_obs, pod_tst = np.zeros(n_pairs, np.int64), np.zeros(n_pairs, np.int64)
    pod_obs[pair_ids] = pod_df.n_observed.values
    pod_tst[pair_ids] = pod_df.n_tested.values
    max_tested = int(pod_tst.max()) + 1
    del pod_df, codes, keys, order

    # ---- load every study once ------------------------------------------------
    with mp.Pool(snakemake.threads) as pool:
        studies = pool.map(load_study, snakemake.input.consituent_studies, chunksize=8)

    # ---- reference fit --------------------------------------------------------
    # only the largest studies have a explicit search space large enough to estimate on
    ref = sorted(studies, key=lambda s: -s[3].sum())[:n_ref]
    ref_obs, ref_tst = np.zeros(n_pairs, np.int64), np.zeros(n_pairs, np.int64)
    for _, ids, obs, tst in ref:
        ref_obs[ids] += obs
        ref_tst[ids] += tst
    m = ref_tst > 0
    fit = fit_mixture(*unique_test_comb(ref_obs[m], ref_tst[m]), max_iter)
    assert not fit["swapped"], "Reference fit swapped components"
    fit["reference_studies"] = ",".join(s[0] for s in ref)

    # ---- external evidence of each study's informative pairs -----------------
    ext = {}
    for name, ids, obs, tst in studies:
        ext_o, ext_t = pod_obs[ids] - obs, pod_tst[ids] - tst
        rep = (obs > 0) & (ext_t > 0)
        ext[name] = ext_o[rep], ext_t[rep]
    ref_cnt = {s[0]: confirm_counts(*ext[s[0]]) for s in ref}
    tot = sum(ref_cnt.values())

    # ---- theta_ref: share of the reference detections that are real -----------
    n_r, t_r = ref_obs[m], ref_tst[m]
    pA = fit["pi"] * betabinom.pmf(n_r, t_r, fit["a"], fit["b"]) # Expected true Positives
    pB = (1 - fit["pi"]) * binom.pmf(n_r, t_r, fit["phi"]) # Expected false Positives
    fit["theta_ref_fit"] = (pA / (pA + pB) * n_r).sum() / n_r.sum()
    # check: the confirmation rate levels off at theta_ref for hits retested >= 10x
    hi = T_BINS >= 5
    q5 = 1 - (1 - fit["phi"]) ** 10 # chance of at least one FP
    print(f"Number of hits retested >= 10x: {tot[1, hi].sum()}")
    print(f"Number of tests retested >= 10x: {tot[0, hi].sum()}")
    fit["n_plateau"] = tot[0, hi].sum()
    fit["theta_ref_plateau"] = (
        (tot[1, hi].sum() / fit["n_plateau"] - q5) / (1 - q5)
        if fit["n_plateau"]
        else np.nan
    )  # nan: no hits retested >= 10x
    if fit["n_plateau"] >= 500: # enough data to trust the plateau estimate
        theta_ref = fit["theta_ref_plateau"]  
    else:
        theta_ref = fit["theta_ref_fit"]
    fit["theta_ref_used"] = theta_ref
    pd.DataFrame([fit]).to_csv(snakemake.output.global_fit, sep="\t", index=False)

    # ---- c(t): confirmation rate of a true hit, corrected for the reference's noise
    q_mid = 1 - (1 - fit["phi"]) ** T_BINS
    conf_rate = lambda N, K: (K / N - (1 - theta_ref) * q_mid) / theta_ref

    # ---- per-study precision, scored against every other study ----------------
    rows = []
    for name, ids, obs, tst in studies:
        N, K = tot - ref_cnt.get(
            name, 0
        )  # leave the scored study out of its own yardstick
        rows.append(
            dict(
                study=name,
                k_s=float(obs.sum()),
                T_s=float(tst.sum()),
                n_reported=int((obs > 0).sum()),
                **score_reported(*ext[name], conf_rate(N, K), fit["phi"]),
            )
        )
    study_df = pd.DataFrame(rows)

    # shrink toward the mean: a study with two overlapping pairs carries no
    # information about its own reliability, and theta at 0 or 1 gives infinite weight
    m_theta = np.nansum(study_df.theta_g_sum) / np.nansum(study_df.n_informative)
    theta = (
        (study_df.theta_g_sum.fillna(0) + k0 * m_theta) / (study_df.n_informative + k0)
    ).clip(1e-6, 1 - 1e-6)

    # p is global; the study's precision enters through the odds ratio
    #   theta/(1-theta) = pi p / ((1-pi) phi)  ->  log(p/phi) = logit(theta) - logit(pi)
    p = fit["p_mean"]
    study_df["theta_s"], study_df["fdr_s"] = theta, 1 - theta
    study_df["w_detected"] = logit(theta) - logit(fit["pi"])
    study_df["phi_s"] = p / np.exp(study_df.w_detected)
    study_df["w_not_detected"] = np.log((1 - p) / (1 - study_df.phi_s))
    study_df.to_csv(snakemake.output.study_fdr, sep="\t", index=False)

    # ---- per-pair posterior ---------------------------------------------------
    log_odds = np.full(n_pairs, logit(fit["pi"]), dtype=np.float64)
    fdr_tst_sum, tst_sum = np.zeros(n_pairs), np.zeros(n_pairs)
    w = study_df.set_index("study")
    for name, ids, obs, tst in studies:
        log_odds[ids] += obs * w.w_detected[name] + (tst - obs) * w.w_not_detected[name]
        fdr_tst_sum[ids] += tst * w.fdr_s[name]
        tst_sum[ids] += tst

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_fdr = fdr_tst_sum / tst_sum

    pd.DataFrame(
        {
            "pair_id": pair_ids,
            "n_tested": pod_tst[pair_ids],
            "n_observed": pod_obs[pair_ids],
            "mean_fdr_per_test": mean_fdr[pair_ids],
            "log_odds": log_odds[pair_ids],
            "posterior": expit(log_odds[pair_ids]),
        }
    ).to_parquet(snakemake.output.pair_posterior, index=False)


if __name__ == "__main__":
    main()
