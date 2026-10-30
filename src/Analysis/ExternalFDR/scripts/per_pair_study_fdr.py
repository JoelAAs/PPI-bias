import os

os.environ["OMP_NUM_THREADS"] = "1"  # one thread per job, the pool is the parallelism

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logit
import multiprocessing as mp


def unique_test_comb(observations, tests):
    k = int(tests.max()).bit_length()
    keys, w = np.unique(
        (observations.astype(np.int64) << k) | tests, return_counts=True
    )
    observed = keys >> k
    tested = keys & ((1 << k) - 1)
    return observed, tested, w.astype(float)

def fit_binomial_mixture(observations, tests, weights, max_iter=10000, tol=1e-8):
    pi, a, b, phi = 0.01, 0.5, 2.5, 4e-4  # prevalence, alpha, beta, FPR
    ll = ll_old = -np.inf
    converged = False
    for i in range(max_iter):
        A = betabinom.pmf(observations, tests, a, b)  # P(n | t, true)
        B = binom.pmf(observations, tests, phi)  # P(n | t, false)
        density = pi * A + (1 - pi) * B  # Mixture density of detection
        true_detection = (
            pi * A
        ) / density  # P(true | n, t)

        w_td = weights * true_detection  # expected true pairs per (n,t)
        w_fd = weights * (1 - true_detection)  # expected false pairs per (n,t)

        pi = w_td.sum() / weights.sum()  # Update prevalence
        phi = (w_fd * observations).sum() / (
            w_fd * tests
        ).sum()  # Update FPR

        lh_n = lambda x: -np.sum(
            w_td
            * betabinom.logpmf(
                observations, tests, *np.exp(x)
            )  # log to keep parameters positive
        )  # Negative log-likelihood of the true component
        a, b = np.exp(
            minimize(lh_n, np.log([a, b]), method="Nelder-Mead").x
        )  # Update detection alpha and beta

        ll = np.sum(weights * np.log(density))
        if ll - ll_old < tol * abs(ll):
            converged = True
            break
        ll_old = ll

    return dict(
        pi=pi,
        a=a,
        b=b,
        phi=phi,
        true_detection_mean=a / (a + b),
        true_detection_var=a * b / ((a + b) ** 2 * (a + b + 1)),
        loglik=ll,
        swapped=a / (a + b) <= phi,  # else the components changed places
        converged=converged,
    )


def pair_keys(ca, cb, k, directed):
    if not directed:
        ca, cb = np.minimum(ca, cb), np.maximum(ca, cb)
    return (ca.astype(np.int64) << k) | cb


def get_pair_indexer(pod_df):
    n_rows = len(pod_df)
    proteins = np.concatenate([pod_df[prot_a].values, pod_df[prot_b].values])
    codes, id_vocab = pd.factorize(proteins, sort=True)
    id_vocab = pd.Index(id_vocab)

    k = (len(id_vocab) - 1).bit_length()
    keys = pair_keys(codes[:n_rows], codes[n_rows:], k, directed)
    order = np.argsort(keys)
    pod_keys, pod_key_pair_ids = keys[order], pod_df.pair_id.values[order]
    assert (
        np.diff(pod_keys) > 0
    ).all(), "Duplicate (bait, prey) in aggregated protein pair table"
    return id_vocab, k, pod_keys, pod_key_pair_ids


def lookup_pair_ids(bait, prey):
    ca, cb = id_vocab.get_indexer(bait), id_vocab.get_indexer(prey)
    key = pair_keys(ca, cb, k, directed)
    pos = np.searchsorted(pod_keys, key).clip(max=len(pod_keys) - 1)
    hit = (ca >= 0) & (cb >= 0) & (pod_keys[pos] == key)
    return np.where(hit, pod_key_pair_ids[pos], -1)


def load_study(path):
    df = pd.read_csv(path, sep="\t", usecols=[prot_a, prot_b, "n_tested", "n_observed"])
    detections = df.n_observed > 0

    detected_df = df[detections]
    multiple_detecting_baits = detected_df.groupby(prot_b)[prot_a].nunique()
    n_detecting_baits = df[prot_b].map(multiple_detecting_baits).fillna(0)
    # If a prey is detected by atleast two baits, this means that reported interaction COULD have been negative data in the same study.
    # Therefore the Hits/Misses isn't so stacked in the direction of interactions
    ev = (n_detecting_baits - detections.astype(int)) >= 1
    df["obs_ev"] = df.n_observed * ev
    df["tst_ev"] = df.n_tested * ev
    df["pair_id"] = lookup_pair_ids(df[prot_a].values, df[prot_b].values)
    df = (
        df[df.pair_id != -1]
        .groupby("pair_id", sort=False)[["n_tested", "n_observed", "tst_ev", "obs_ev"]]
        .sum()
    )
    return (
        os.path.basename(path).removesuffix(".csv"),
        df.index.values.astype(np.int32),
        df.n_observed.values.astype(np.int32),
        df.n_tested.values.astype(np.int32),
        df.obs_ev.values.astype(np.int32),
        df.tst_ev.values.astype(np.int32),
    )
    


def score_reported(obs, tests, c, phi, tested_bins):
    if len(obs) == 0: # no overlapping, external evidance
        return dict(theta_g_sum=np.nan, n_informative=0, mean_t_ext=np.nan)
    
    conf = obs > 0
    q = 1 - (1 - phi) ** tests # chance of at least one FP
    ct = np.clip(c[np.digitize(tests, tested_bins) - 1], q + 1e-12, 1 - 1e-12)
    likelihood = lambda theata: np.where(
        conf,
        theata * ct + (1 - theata) * q, # chance of observing external confirmation
        1 - theata * ct - (1 - theata) * q) # chance of not confimed
    theta = minimize_scalar(
        lambda th: -np.log(likelihood(th)).sum(), bounds=(0, 1), method="bounded"
    ).x
    g = theta * np.where(conf, ct, 1 - ct) / likelihood(theta)  # P(true | external evidence) = P(Evidance|true)/P(evidance)
    return dict(
        theta_g_sum=g.sum(),
        n_informative=len(obs),
        mean_t_ext=tests.mean(),
        f_confirmed=conf.mean(),
        f_expected_true=ct.mean(),
    )

def confirmed_interaction_rate(tested, detected, summed_q, theta_ref):
    mean_q = summed_q / tested
    c = (detected / tested - (1 - theta_ref) * mean_q) / theta_ref
    return c


def main():
    global directed, prot_a, prot_b
    global id_vocab, k, pod_keys, pod_key_pair_ids

    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    max_iter = snakemake.params.n_em_iterations
    k0 = snakemake.params.shrinkage_pseudocount
    n_ref = snakemake.params.n_reference_studies
    n_bins = 4

    directed = {"undirectional": False, "directional": True}[
        snakemake.wildcards.network_type
    ]
    pod_df = pd.read_parquet(
        snakemake.input.pod_file,
        columns=[prot_a, prot_b, "pair_id", "n_tested", "n_observed"],
    )

    # Indexer
    id_vocab, k, pod_keys, pod_key_pair_ids = get_pair_indexer(pod_df)

    # Aggregate for final edge output
    pair_ids = pod_df.pair_id.values
    n_pairs = int(pair_ids.max()) + 1  # arrays are indexed by pair_id
    assert (
        np.bincount(pair_ids, minlength=n_pairs).max() <= 1
    ), "Duplicate pair_id in POD"
    pod_obs = np.zeros(n_pairs, np.int64)
    pod_tst = np.zeros(n_pairs, np.int64)
    pod_obs[pair_ids] = pod_df.n_observed.values
    pod_tst[pair_ids] = pod_df.n_tested.values
    del pod_df
    # Load all studies
    with mp.Pool(snakemake.threads) as pool:
        studies = pool.map(load_study, snakemake.input.consituent_studies, chunksize=8)

    # Define reference set
    ref = sorted(studies, key=lambda s: -s[3].sum())[
        :n_ref
    ]  # the N studies with the largest amount of tests
    ref_obs = np.zeros(n_pairs, np.int64)
    ref_tst = np.zeros(n_pairs, np.int64)
    for _, ids, _, _, obs_ev, tst_ev in ref:
        ref_obs[ids] += obs_ev
        ref_tst[ids] += tst_ev
    ref_pair_ids = ref_tst > 0
    fit = fit_binomial_mixture(
        *unique_test_comb(ref_obs[ref_pair_ids], ref_tst[ref_pair_ids]),
        max_iter)
    assert not fit["swapped"], "Reference fit swapped components"
    fit["reference_studies"] = ";".join(s[0] for s in ref)

    # Define probability of confirming a True interaction within N tests
    tested_bins = np.arange(1, n_bins + 1)  # bins of external tests

    def _confirm_counts(observations, tests):
        pair_tested_bins = np.digitize(tests, tested_bins) - 1
        return np.stack(
            [
                np.bincount(pair_tested_bins, minlength=len(tested_bins)),
                np.bincount(
                    pair_tested_bins, weights=observations > 0, minlength=len(tested_bins)),
                np.bincount(
                    pair_tested_bins, weights=1 - (1 - fit["phi"]) ** tests, minlength=len(tested_bins)),

            ]
        )
    
    ref_names = {s[0] for s in ref}
    ext = {}
    for name, ids, obs, tst, obs_ev, tst_ev in studies:
        own = name in ref_names  # a reference study does not confirm itself
        ext_o = ref_obs[ids] - (obs_ev if own else 0)
        ext_t = ref_tst[ids] - (tst_ev if own else 0)
        rep = (obs > 0) & (ext_t > 0)  # the study's own hits that a reference study also tested
        ext[name] = ext_o[rep], ext_t[rep]
        
    ref_cnt = {
            s[0]: _confirm_counts(*ext[s[0]]) for s in ref
        }
    tot = sum(ref_cnt.values()) 
    
    # theta_ref: share of the reference detections that come from true pairs (mixture fit)
    n_r = ref_obs[ref_pair_ids]
    t_r = ref_tst[ref_pair_ids]
    pA = fit["pi"] * betabinom.pmf(n_r, t_r, fit["a"], fit["b"])  # prior x P(counts | true)
    pB = (1 - fit["pi"]) * binom.pmf(n_r, t_r, fit["phi"])  # prior x P(counts | false)
    fit["theta_ref_fit"] = (pA / (pA + pB) * n_r).sum() / n_r.sum()
    theta_ref = fit["theta_ref_fit"]

    # Output global fit data
    pd.DataFrame([fit]).to_csv(snakemake.output.global_fit, sep="\t", index=False)

    # Score studies
    rows = []
    for name, ids, obs, tst, *_ in studies:
        external_tested_seen = tot - ref_cnt.get(name, 0)  # leave a reference study out
        rows.append(
            dict(
                study=name,
                k_s=float(obs.sum()),
                T_s=float(tst.sum()),
                n_reported=int((obs > 0).sum()),
                **score_reported(
                    *ext[name],
                    confirmed_interaction_rate(*external_tested_seen, theta_ref),
                    fit["phi"],
                    tested_bins),
            )
        )
        
    study_df = pd.DataFrame(rows)

    m_theta = np.nansum(study_df.theta_g_sum) / np.nansum(study_df.n_informative)
    theta = (
        (study_df.theta_g_sum.fillna(0) + k0 * m_theta) / (study_df.n_informative + k0)
    ).clip(1e-6, 1 - 1e-6)
    
    td_mean = fit["true_detection_mean"]
    study_df["theta_s"]    = theta
    study_df["fdr_s"]      = 1 - theta
    study_df["w_detected"] = logit(theta) - logit(fit["pi"])
    study_df["phi_s"]      = td_mean / np.exp(study_df.w_detected)
    study_df["w_not_detected"] = np.log((1 - td_mean) / (1 - study_df.phi_s))
    study_df.to_csv(snakemake.output.study_fdr, sep="\t", index=False)
    
    log_odds = np.full(n_pairs, logit(fit["pi"]), dtype=np.float64)
    fdr_tst_sum, tst_sum = np.zeros(n_pairs), np.zeros(n_pairs)
    w = study_df.set_index("study")
    for name, ids, obs, tst, *_ in studies:  # full counts: every hit counts for its pair
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