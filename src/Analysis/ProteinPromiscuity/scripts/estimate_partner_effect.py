from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats


def _eval_rho(r, role_pairs, role, rng, cutoff, n_permutations):
    p_val_distribution = role_pairs.groupby(role).apply(uniformity_p, r, n_permutations, rng, include_groups=False)

    p_val_distribution = p_val_distribution[p_val_distribution > cutoff] # as there might be those with rejected H0
    p_val_distribution = p_val_distribution[~p_val_distribution.isna()]
    #rescale to 0-1
    p_val_distribution = (p_val_distribution - cutoff) / (1 - cutoff)

    if len(p_val_distribution) < 20:
        # To few points are bad.
        W,p = np.inf, 1
    else:
        cvm = stats.cramervonmises(p_val_distribution, 'uniform')
        W, p = cvm.statistic, cvm.pvalue
    return {"W": W, "r0": r, "p": p, "n": len(p_val_distribution)}


def estimate_rho(unknown_role_pairs, role, tot, rng, threads, n_permutations):
    r0 = 0.3
    max_iter = 100
    cutoff = 0.3
    max_iter = 3
    steps = 10
    role_pairs = unknown_role_pairs[unknown_role_pairs[role].isin(tot.index[tot["D"].notna()])]

    def _eval_linspace(rho):
        tasks = [(r, role_pairs, role, c_rng, cutoff, n_permutations) for r, c_rng in zip(rho, rng.spawn(len(rho)))]
        with Pool(threads) as pool:
            results = pool.starmap(_eval_rho, tasks)
        best = {"W": np.inf, "r0": None, "p": 0.0, "n": 0}
        for res in results:
            if res["W"] < best["W"]:
                best = res
        return best

    rho = np.linspace(0.01, 0.99, steps)
    for i in range(max_iter):
        best = _eval_linspace(rho)
        
        if best["p"] > 0.05:
            return best["r0"]
        else:
            r0_ = best["r0"]
            idx = np.searchsorted(rho, r0_)
            max_idx = min(idx + 1, len(rho) - 1)
            min_idx = max(idx - 1, 0)
            rho = np.linspace(rho[min_idx], rho[max_idx], steps)
    return best["r0"]    

def uniformity_p(g, r, n_permutations, rng):
    n = g["n"].to_numpy()
    p = g["p"].iloc[0]
    B = n_permutations
    
    sim = rng.binomial(n, rng.beta(p * (1 - r) / r, (1 - p) * (1 - r) / r, size=(B, n.size)))
    q = sim.sum(axis=1, keepdims=True) / n.sum()  # re-estimate p as the observed stat does
    chi = ((sim - n * q) ** 2 / np.maximum(n * q * (1 - q), 1e-12)).sum(axis=1)
    return (1 + (chi <= g["chi"].sum()).sum()) / (B + 1)


def main():
    n_permutations = snakemake.params.n_simulations
    role = snakemake.wildcards.role
    partner = "prey" if role == "bait" else "bait"
    df = pd.read_csv(snakemake.input.row_wise_data, sep="\t")
    pair = df.groupby([role, partner], as_index=False).agg(
        n=("detection", "size"), k=("detection", "sum")
    )

    pair = pair[pair["n"] >= 2].copy()     # n==1 pairs are uninformative by construction
    pair["det"] = pair["k"] > 0            # detected at least once
    pair["full"] = pair["k"] == pair["n"]  # detected on every test
    tot = pair.groupby(role).agg(
        n_pairs=("n", "size"), n_obs=("n", "sum"), n_pos=("k", "sum"),
        n_det=("det", "sum"), n_full=("full", "sum"),
    )
    tot["p"] = tot["n_pos"] / tot["n_obs"]
    pair = pair.join(tot["p"], on=role)
    e = pair["n"] * pair["p"]
    v = pair["n"] * pair["p"] * (1 - pair["p"])
    pair["chi"] = ((pair["k"] - e) ** 2 / v).where(v > 0, 0.0)
    tot["D"] = pair.groupby(role)["chi"].sum() / (tot["n_pairs"] - 1)
    tot["f"] = (tot["n_full"] / tot["n_det"]).where(tot["n_det"] > 0) # Fraction of full replication

    # p = 0 or 1: every pair is uniform by construction, so D carries no information
    tot.loc[(tot["p"] == 0) | (tot["p"] == 1), "D"] = np.nan
    tot = tot[tot["n_pairs"] >= snakemake.params.min_pairs]

    rng = np.random.default_rng(0)

    best_rho = estimate_rho(
        pair,
        role,
        tot,
        rng,
        threads=snakemake.threads,
        n_permutations=int(n_permutations/2))
    
    tot["p_vs_specific"] = (
        pair[pair[role].isin(tot.index[tot["D"].notna()])].groupby(role).apply(uniformity_p, best_rho, n_permutations, rng, include_groups=False)
    )
    tot["selected_rho"] = best_rho
    detected = tot["n_pos"] > 0
    tot["score"] = np.nan
    tot.loc[detected, "score"] = (-tot.loc[detected, "D"]).rank(
        pct=True
    ) + tot.loc[detected, "p"].rank(pct=True)
    tot.sort_values("score", ascending=False, na_position="last").to_csv(
        snakemake.output.ranked, sep="\t"
    )

if __name__ == "__main__":
    main()