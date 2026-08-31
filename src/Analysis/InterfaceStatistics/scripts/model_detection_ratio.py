"""Detection rate vs. interface size, on pairs with a matched PDB interface.

    n_observed ~ Binomial(n_tested, p)
    logit(p) = b0 + b1*log(interface_residues) + b2*log(min_length)
                  + b3*log(n_tested)

b1 is the estimand. b2 removes the length -> interface-size and
length -> MS-detectability path; b3 removes the study-attention path.

Dependence between pairs sharing a protein, and overdispersion between
tests of the same pair, are handled by bootstrapping over proteins rather
than by random effects -- the CI comes from the bootstrap, not the GLM.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

TERMS = ["log_interface", "log_min_length", "log_n_tested"]


def read_fasta_lengths(fasta_file):
    lengths, acc, n = {}, None, 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if acc:
                    lengths[acc] = n
                acc, n = line.split("|")[1], 0
            else:
                n += len(line.strip())
    if acc:
        lengths[acc] = n
    return lengths


def prepare(df, lengths):
    df = df.copy()
    df["protein_a"] = df[["uniprot_id_bait", "uniprot_id_prey"]].min(axis=1)
    df["protein_b"] = df[["uniprot_id_bait", "uniprot_id_prey"]].max(axis=1)

    df["min_length"] = np.minimum(
        df["protein_a"].map(lengths), df["protein_b"].map(lengths)
    )
    df = df.dropna(subset=["min_length"])
    df = df[(df.interface_residues > 0) & (df.n_tested > 0)]
    for col, src in [("log_interface", "interface_residues"),
                     ("log_min_length", "min_length"),
                     ("log_n_tested", "n_tested")]:
        df[col] = np.log(df[src]) - np.log(df[src]).mean()
    return df.reset_index(drop=True)


def fit(df):
    X = sm.add_constant(df[TERMS])
    y = np.column_stack([df.n_observed, df.n_tested - df.n_observed])
    return sm.GLM(y, X, family=sm.families.Binomial()).fit()


# def bootstrap_b1(df, n_boot, rng):
#     """Resample proteins, not tests: a pair is drawn when either of its
#     proteins is, which propagates the shared-protein dependence."""
#     proteins = pd.unique(pd.concat([df.protein_a, df.protein_b]))
#     rows = {p: np.flatnonzero((df.protein_a == p) | (df.protein_b == p))
#             for p in proteins}
#     out = []
#     for _ in range(n_boot):
#         drawn = rng.choice(proteins, size=len(proteins), replace=True)
#         idx = np.concatenate([rows[p] for p in drawn])
#         try:
#             out.append(fit(df.iloc[idx].reset_index(drop=True)).params["log_interface"])
#         except Exception:
#             pass
#     return np.array(out)

def bootstrap_b1(df, n_boot, rng):
    """Resample proteins and refit the model, returning all coefficients."""
    proteins = pd.unique(pd.concat([df.protein_a, df.protein_b]))
    rows = {
        p: np.flatnonzero(
            (df.protein_a == p) | (df.protein_b == p)
        )
        for p in proteins
    }

    out = []
    for _ in range(n_boot):
        drawn = rng.choice(
            proteins,
            size=len(proteins),
            replace=True
        )
        idx = np.concatenate([rows[p] for p in drawn])

        try:
            params = fit(
                df.iloc[idx].reset_index(drop=True)
            ).params
            out.append(params)
        except Exception:
            pass

    return pd.DataFrame(out)

def permute_b1(df, n_perm, rng):
    """Shuffle interface size within min_length deciles, so the null keeps
    the length relationship and only breaks the interface signal."""
    decile = pd.qcut(df.log_min_length, 10, labels=False, duplicates="drop")
    out = []
    for _ in range(n_perm):
        d = df.copy()
        d["log_interface"] = (
            d.groupby(decile)["log_interface"]
            .transform(lambda x: rng.permutation(x.values))
        )
        out.append(fit(d).params["log_interface"])
    return np.array(out)


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    n_boot = snakemake.params.get("n_bootstrap", 500)
    n_perm = snakemake.params.get("n_permutations", 200)

    log = open(snakemake.log[0], "w")

    lengths = read_fasta_lengths(snakemake.input.gene_fasta)
    df = prepare(pd.read_parquet(snakemake.input.interface_annotated), lengths)
    print(f"{len(df)} pairs, {df.n_tested.sum()} tests, "
          f"{df.n_observed.sum()} observations", file=log, flush=True)

    m = fit(df)
    b1 = m.params["log_interface"]
    boot = bootstrap_b1(df, n_boot, rng)
    perm = permute_b1(df, n_perm, rng)

    to_or = lambda b: np.exp(b * np.log(2))
    or_est = to_or(b1)
    or_lo, or_hi = to_or(np.percentile(boot, [2.5, 97.5]))
    perm_p = (np.sum(np.abs(perm) >= abs(b1)) + 1) / (len(perm) + 1)

    coefs = pd.DataFrame({
        "term": m.params.index,
        "estimate": m.params.values,
        "se_naive": m.bse.values,
        "p_value_naive": m.pvalues.values,
        "or_per_doubling": to_or(m.params.values),
    })

    # Add bootstrap percentile CIs
    boot_lo = []
    boot_hi = []

    for term in m.params.index:
        vals = boot[term].dropna()
        lo, hi = np.percentile(vals, [2.5, 97.5])

        boot_lo.append(to_or(lo))
        boot_hi.append(to_or(hi))

    coefs["or_bootstrap_low"] = boot_lo
    coefs["or_bootstrap_high"] = boot_hi
    coefs.to_csv(snakemake.output.coefficients, sep="\t", index=False)

    print(m.summary(), file=log)
    print(f"\nOR per doubling of interface size: {or_est:.3f} "
          f"[{or_lo:.3f}, {or_hi:.3f}] (protein bootstrap, n={len(boot["log_interface"])})",
          file=log)
    print(f"Naive GLM SE understates the bootstrap SD by "
          f"{boot["log_interface"].std() / m.bse['log_interface']:.2f}x", file=log)
    print(f"Permutation null b1: mean {perm.mean():.4f}, sd {perm.std():.4f}; "
          f"empirical two-sided p = {perm_p:.4f}", file=log)

    # partial effect, other covariates held at their geometric means
    grid = np.linspace(df.log_interface.min(), df.log_interface.max(), 200)
    Xg = sm.add_constant(pd.DataFrame({
        "log_interface": grid, "log_min_length": 0.0, "log_n_tested": 0.0
    }))
    pred = m.get_prediction(Xg).summary_frame()
    x = np.exp(grid + np.log(df.interface_residues).mean())

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.fill_between(x, pred["mean_ci_lower"], pred["mean_ci_upper"],
                    color="firebrick", alpha=0.2)
    ax.plot(x, pred["mean"], color="firebrick")
    ax.set_xscale("log")
    ax.set_xlabel("Interface size (residues)")
    ax.set_ylabel("Fitted detection probability")
    ax.set_title(f"OR {or_est:.2f} per doubling [{or_lo:.2f}-{or_hi:.2f}]")
    fig.tight_layout()
    fig.savefig(snakemake.output.partial_effect, dpi=200)

    np.savez(snakemake.output.resamples, bootstrap=boot, permutation=perm)
    print("Done.", file=log)
    log.close()