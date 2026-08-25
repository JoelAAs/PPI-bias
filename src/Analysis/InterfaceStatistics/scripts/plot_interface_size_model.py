import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta


if __name__ == "__main__":
    model_file = snakemake.input.model
    protein_lengths_file = snakemake.input.protein_lengths
    output_file = snakemake.output.plot
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    bundle = joblib.load(model_file)
    m, phi = bundle["model"], bundle["phi"]

    df = pd.read_csv(protein_lengths_file)
    print(f"Loaded {len(df)} fitted pairs", file=log, flush=True)

    L_grid = np.logspace(np.log10(df.min_length.min()), np.log10(df.min_length.max()), 200)
    mu = m.predict(pd.DataFrame({"logL": np.log(L_grid)}))
    a, b = mu * phi, (1 - mu) * phi
    lower = beta.ppf(0.025, a, b)
    upper = beta.ppf(0.975, a, b)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df.min_length, df.r, s=8, alpha=0.3, color="grey", label="observed pairs")
    ax.plot(L_grid, mu, color="C0", label="fitted mean")
    ax.fill_between(L_grid, lower, upper, color="C0", alpha=0.2, label="95% credible interval")
    ax.set_xscale("log")
    ax.set_xlabel("min(length_bait, length_prey)")
    ax.set_ylabel("interface_residues / min_length (r)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_file, dpi=200)

    print("Done.", file=log, flush=True)
    log.close()
