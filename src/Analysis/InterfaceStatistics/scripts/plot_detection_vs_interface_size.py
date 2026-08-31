import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


if __name__ == "__main__":
    interface_annotated_file = snakemake.input.interface_annotated
    output_file = snakemake.output.plot
    log_file = snakemake.log[0]
    n_bins = 12

    log = open(log_file, "w")

    df = pd.read_parquet(interface_annotated_file)
    df = df[df["n_tested"] > 0].copy()
    df["detection_ratio"] = df["n_observed"] / df["n_tested"]
    print(f"Loaded {len(df)} pairs with interface size and detection counts", file=log, flush=True)

    # pooled (sum n_observed / sum n_tested) rather than mean-of-ratios,
    # so bins aren't dominated by many low-n_tested pairs stuck at 0 or 1
    edges = np.logspace(np.log10(df.interface_residues.min()), np.log10(df.interface_residues.max()), n_bins + 1)
    df["size_bin"] = pd.cut(df.interface_residues, edges, include_lowest=True)
    binned = df.groupby("size_bin", observed=True).agg(
        n_observed=("n_observed", "sum"),
        n_tested=("n_tested", "sum"),
        mid=("interface_residues", "median"),
    )
    binned["pooled_ratio"] = binned["n_observed"] / binned["n_tested"]
    binned = binned.dropna(subset=["pooled_ratio"])
    print(f"Binned into {len(binned)} interface-size bins", file=log, flush=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        df.interface_residues, df.detection_ratio,
        c=df.n_tested, norm=LogNorm(), cmap="viridis",
        s=10, alpha=0.4, edgecolors="none", label="pairs",
    )
    ax.plot(binned.mid, binned.pooled_ratio, color="red", marker="o", linewidth=2, label="pooled ratio per bin")
    ax.set_xscale("log")
    ax.set_xlabel("Interface size (residues, log scale)")
    ax.set_ylabel("Detection ratio (n_observed / n_tested)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("n_tested")
    fig.tight_layout()
    fig.savefig(output_file, dpi=200)

    print("Done.", file=log, flush=True)
    log.close()
