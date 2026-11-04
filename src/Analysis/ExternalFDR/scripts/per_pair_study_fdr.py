import os

os.environ["OMP_NUM_THREADS"] = "1"  # one thread per job, the pool is the parallelism

# THIS WILL ONLY WORK FOR MS


import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logit
import multiprocessing as mp


def fit_binomial_mixture(observations, tests, weights, max_iter=10000, tol=1e-8, init=(0.01, 0.5, 2.5, 4e-4)):
    pi, a, b, phi = init  # prevalence, alpha, beta, FPR
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


def pair_keys(ca, cb, k):

    return (ca.astype(np.int64) << k) | cb


def get_pair_indexer(pod_df):
    n_rows = len(pod_df)
    proteins = np.concatenate([pod_df[prot_a].values, pod_df[prot_b].values])
    codes, id_vocab = pd.factorize(proteins, sort=True)
    id_vocab = pd.Index(id_vocab)

    k = (len(id_vocab) - 1).bit_length()
    keys = pair_keys(codes[:n_rows], codes[n_rows:], k)
    order = np.argsort(keys)
    pod_keys, pod_key_pair_ids = keys[order], pod_df.pair_id.values[order]
    assert (
        np.diff(pod_keys) > 0
    ).all(), "Duplicate (bait, prey) in aggregated protein pair table"
    return id_vocab, k, pod_keys, pod_key_pair_ids


def lookup_pair_ids(bait, prey):
    ca, cb = id_vocab.get_indexer(bait), id_vocab.get_indexer(prey)
    key = pair_keys(ca, cb, k)
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


def padded_reference(ref, pair_bait, pair_prey, n_proteins):
    """
    Reference studies on one fixed prey universe U, the preys in the search space of the largest reference study.
    Every bait of a reference study counts as tested against every prey in U (self pairs excluded), so whether a
    test counts never depends on its outcome and all reference tests share one prior pi.
    :param ref: reference studies as returned by load_study, largest first
    :param pair_bait: (np.ndarray) bait protein code per pair_id
    :param pair_prey: (np.ndarray) prey protein code per pair_id
    :param n_proteins: (int) size of the protein vocabulary
    :return: (dict) pair_bait, pair_prey as given
                    in_U: per protein, True if in the prey universe
                    n_ref_baits: per protein, number of reference studies that used it as bait
                    ref_det: per pair_id, number of reference studies detecting it (prey in U, no self pairs)
    """
    in_U = np.zeros(n_proteins, bool)
    in_U[pair_prey[ref[0][1]]] = True
    n_ref_baits = np.zeros(n_proteins, np.int64)
    ref_det = np.zeros(len(pair_bait), np.int8)
    for _, ids, obs, _ in ref:
        n_ref_baits[np.unique(pair_bait[ids])] += 1
        det = ids[obs > 0]
        ref_det[det[in_U[pair_prey[det]] & (pair_bait[det] != pair_prey[det])]] += 1
    return dict(pair_bait=pair_bait, pair_prey=pair_prey, in_U=in_U, n_ref_baits=n_ref_baits, ref_det=ref_det)


def reference_histogram(rs):
    """
    Padded reference pairs counted per (n_observed, n_tested), the input of the mixture fit. Never detected
    pairs are counted per bait, not listed
    :param rs: (dict) from padded_reference
    :return: (observed, tested, weights) arrays, one entry per occurring (n_observed, n_tested)
    """
    in_U, n_ref_baits = rs["in_U"], rs["n_ref_baits"]
    pairs_per_bait = in_U.sum() - in_U  # self pair excluded
    det_ids = np.flatnonzero(rs["ref_det"])
    det_n = rs["ref_det"][det_ids].astype(np.int64)
    det_t = n_ref_baits[rs["pair_bait"][det_ids]]
    t_max = int(n_ref_baits.max())
    counts = np.zeros((t_max + 1, t_max + 1))  # [n_observed, n_tested]
    np.add.at(counts, (det_n, det_t), 1)
    counts[0] = np.bincount(n_ref_baits, weights=pairs_per_bait, minlength=t_max + 1) - counts[1:].sum(0)
    counts[0, 0] = 0  # proteins no reference study used as bait
    observed, tested = np.nonzero(counts)
    return observed, tested, counts[observed, tested]


def external_counts(ids, obs, own, rs):
    """
    Padded reference evidence for the pairs a study tested: a reference study tested a pair if it used its bait
    and the prey is in U
    :param ids: (np.ndarray) pair_ids the study tested
    :param obs: (np.ndarray) the study's detections of those pairs
    :param own: (bool) the study is itself a reference study, left out of its own evidence
    :param rs: (dict) from padded_reference
    :return: (ext_o, ext_t) number of reference studies detecting / testing each pair
    """
    b, p = rs["pair_bait"][ids], rs["pair_prey"][ids]
    padded = rs["in_U"][p] & (b != p)
    ext_t = rs["n_ref_baits"][b] * padded - (own & padded)
    ext_o = rs["ref_det"][ids].astype(np.int64) - (own & padded & (obs > 0))
    return ext_o, ext_t


def true_posterior(fit, n, t):
    """P(true | n of t reference tests detected) under the mixture fit"""
    pA = fit["pi"] * betabinom.pmf(n, t, fit["a"], fit["b"])  # prior x P(counts | true)
    pB = (1 - fit["pi"]) * binom.pmf(n, t, fit["phi"])  # prior x P(counts | false)
    return pA / (pA + pB)


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


def diagnose_reference_fit(fit, reference, studies, ref_names, rs, study_df, max_iter):
    """
    Checks behind the theta_s >= pi floor, written to the three diagnostics outputs:
      1. refits: is pi stable when the reference mixture is refitted from other starting values
      2. predictive_check: does the fitted mixture reproduce the reference (n_observed | n_tested) counts
      3. floor_sensitivity: how many detections sit in studies the floor silences, near the floor,
         or on the wrong side of their own baseline pi_s
    :param reference: (observed, tested, weights) from reference_histogram
    :param rs: (dict) from padded_reference
    :return: (np.ndarray) pi_s per study (same order as studies), the mean P(true | external reference counts)
             over all pairs the study tested that the padded reference also tested; the study's own counts left out
    """
    pi, a, b, phi = fit["pi"], fit["a"], fit["b"], fit["phi"]
    n_u, t_u, w_u = reference

    # 1. Refits, first row is the default start and should reproduce the main fit
    starts = [
        (0.01, 0.5, 2.5, 4e-4),  # default
        (0.001, 0.5, 2.5, 4e-4),  # rarer interactions
        (0.1, 0.5, 2.5, 4e-4),  # commoner interactions
        (0.01, 1.0, 1.0, 4e-4),  # uniform detectability, mean 0.5
        (0.01, 5.0, 45.0, 4e-4),  # homogeneous low detectability, mean 0.1
        (0.01, 0.5, 2.5, 4e-3),  # noisier false detections
    ]
    with mp.Pool(min(snakemake.threads, len(starts))) as pool:
        refits = pool.starmap(
            fit_binomial_mixture, [(n_u, t_u, w_u, max_iter, 1e-8, s) for s in starts]
        )
    pd.DataFrame(
        [dict(zip(["start_pi", "start_a", "start_b", "start_phi"], s), **r) for s, r in zip(starts, refits)]
    ).to_csv(snakemake.output.refits, sep="\t", index=False)

    # 2. Observed vs expected number of reference pairs per (n_tested, n_observed)
    ts = np.unique(t_u)
    t_rep = np.repeat(ts, ts + 1)
    n_rep = np.concatenate([np.arange(t + 1) for t in ts])
    n_pairs_t = pd.Series(w_u).groupby(t_u).sum().reindex(t_rep).values
    observed = (
        pd.Series(w_u, index=pd.MultiIndex.from_arrays([t_u, n_u]))
        .reindex(pd.MultiIndex.from_arrays([t_rep, n_rep]), fill_value=0)
        .values
    )
    ppc = pd.DataFrame(
        dict(
            n_tested=t_rep,
            n_observed=n_rep,
            observed=observed,
            expected_true=n_pairs_t * pi * betabinom.pmf(n_rep, t_rep, a, b),
            expected_false=n_pairs_t * (1 - pi) * binom.pmf(n_rep, t_rep, phi),
        )
    )
    ppc["expected"] = ppc.expected_true + ppc.expected_false
    ppc.to_csv(snakemake.output.predictive_check, sep="\t", index=False)

    # 3a. pi_s: the baseline each study's hits have to beat
    T = int(t_u.max())
    t_g, n_g = np.tril_indices(T + 1)  # every n <= t
    post = np.full((T + 1, T + 1), np.nan)
    post[n_g, t_g] = true_posterior(fit, n_g, t_g)

    pi_s = np.full(len(studies), np.nan)
    for i, (name, ids, obs, tst) in enumerate(studies):
        ext_o, ext_t = external_counts(ids, obs, name in ref_names, rs)  # a reference study does not set its own baseline
        tested = ext_t > 0
        if tested.any():  # no overlap with the reference, no baseline (as in score_reported)
            pi_s[i] = post[ext_o[tested], ext_t[tested]].mean()

    # 3b. Share of detections affected by where the floor sits
    theta, k = study_df.theta_s.values, study_df.k_s.values
    groups = {
        "all": np.ones(len(theta), bool),
        "below_pi": theta < pi,  # silenced by the floor
        "near_floor": (theta >= pi / 2) & (theta <= 2 * pi),  # would flip with a 2x different pi
        "silenced_above_pi_s": (theta > pi_s) & (theta < pi),  # beat their own baseline, silenced anyway
        "credited_below_pi_s": (theta > pi) & (theta < pi_s),  # lost to their own baseline, credited anyway
    }
    pd.DataFrame(
        [
            dict(group=g, n_studies=int(m.sum()), n_detections=k[m].sum(), share_detections=k[m].sum() / k.sum())
            for g, m in groups.items()
        ]
    ).to_csv(snakemake.output.floor_sensitivity, sep="\t", index=False)

    return pi_s


def main():
    global prot_a, prot_b
    global id_vocab, k, pod_keys, pod_key_pair_ids

    prot_a = f"{snakemake.params.id_pattern}_bait"
    prot_b = f"{snakemake.params.id_pattern}_prey"
    max_iter = snakemake.params.n_em_iterations
    k0 = snakemake.params.shrinkage_pseudocount
    n_ref = snakemake.params.n_reference_studies
    n_bins = 4

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
    pair_bait = np.full(n_pairs, -1, np.int32)  # protein codes per pair_id
    pair_prey = np.full(n_pairs, -1, np.int32)
    pair_bait[pod_key_pair_ids] = pod_keys >> k
    pair_prey[pod_key_pair_ids] = pod_keys & ((1 << k) - 1)

    # Load all studies
    with mp.Pool(snakemake.threads) as pool:
        studies = pool.map(load_study, snakemake.input.consituent_studies, chunksize=8)

    # Define reference set
    ref = sorted(studies, key=lambda s: -s[3].sum())[
        :n_ref
    ]  # the N studies with the largest amount of tests
    rs = padded_reference(ref, pair_bait, pair_prey, len(id_vocab))  # prey universe: the largest study's preys
    reference = reference_histogram(rs)
    fit = fit_binomial_mixture(*reference, max_iter)
    assert not fit["swapped"], "Reference fit swapped components"
    fit["reference_studies"] = ";".join(s[0] for s in ref)
    fit["prey_universe"] = ref[0][0]
    fit["n_prey_universe"] = int(rs["in_U"].sum())

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
    for name, ids, obs, tst in studies:
        ext_o, ext_t = external_counts(ids, obs, name in ref_names, rs)  # a reference study does not confirm itself
        rep = (obs > 0) & (ext_t > 0)  # the study's own hits that a reference study also tested
        ext[name] = ext_o[rep], ext_t[rep]
        
    ref_cnt = {
            s[0]: _confirm_counts(*ext[s[0]]) for s in ref
        }
    tot = sum(ref_cnt.values()) 
    
    # theta_ref: share of the reference detections that come from true pairs (mixture fit)
    n_r, t_r, w_r = reference
    fit["theta_ref_fit"] = (w_r * true_posterior(fit, n_r, t_r) * n_r).sum() / (w_r * n_r).sum()
    theta_ref = fit["theta_ref_fit"]

    # Output global fit data
    pd.DataFrame([fit]).to_csv(snakemake.output.global_fit, sep="\t", index=False)

    # Score studies
    rows = []
    for name, ids, obs, tst in studies:
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
    study_df["w_detected"] = np.maximum(logit(theta) - logit(fit["pi"]), 0) # if a study is worse than random, its uninformative, not, anti-evidence
    study_df["phi_s"]      = td_mean / np.exp(study_df.w_detected)
    study_df["w_not_detected"] = np.log((1 - td_mean) / (1 - study_df.phi_s))
    study_df["pi_s"] = diagnose_reference_fit(
        fit, reference, studies, ref_names, rs, study_df, max_iter
    )
    study_df.to_csv(snakemake.output.study_fdr, sep="\t", index=False)

    llr_det, llr_nd = np.zeros(n_pairs), np.zeros(n_pairs)  # summed log Bayes factors of detections / misses
    fdr_tst_sum, tst_sum = np.zeros(n_pairs), np.zeros(n_pairs)
    w = study_df.set_index("study")
    for name, ids, obs, tst in studies:  # full counts: every hit counts for its pair
        llr_det[ids] += obs * w.w_detected[name]
        llr_nd[ids] += (tst - obs) * w.w_not_detected[name]
        fdr_tst_sum[ids] += tst * w.fdr_s[name]
        tst_sum[ids] += tst
    log_odds = logit(fit["pi"]) + llr_det + llr_nd

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_fdr = fdr_tst_sum / tst_sum

    pair_df = pd.read_parquet(
        snakemake.input.pod_file, columns=[prot_a, prot_b, "pubmed_id"]
    )
    pair_df["n_observed"] = pod_obs[pair_ids]
    pair_df["n_tested"] = pod_tst[pair_ids]
    pair_df["mean_fdr_per_test"] = mean_fdr[pair_ids]
    pair_df["llr_detected"] = llr_det[pair_ids]
    pair_df["llr_not_detected"] = llr_nd[pair_ids]
    pair_df["log_odds"] = log_odds[pair_ids]
    pair_df["posterior"] = expit(log_odds[pair_ids])
    pair_df.to_parquet(snakemake.output.pair_posterior, index=False)


if __name__ == "__main__":
    main()