import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


def bracket(ax, x1, x2, y, pval):
    label = "ns" if pval >= 0.05 else ("*" if pval >= 0.01 else ("**" if pval >= 0.001 else "***"))
    ax.plot([x1, x1, x2, x2], [y, y * 1.02, y * 1.02, y], lw=1.2, c="black")
    ax.text((x1 + x2) / 2, y * 1.02, label, ha="center", va="bottom")


def main():
    pe = pd.read_csv(snakemake.input.partner_effect, sep="\t").rename(
        columns={snakemake.wildcards.role: "protein_id"}).dropna(subset=["p_vs_specific"])
    selected = pd.read_csv(snakemake.input.selected_proteins, sep="\t").rename(
        columns={"UniProt": "uniprot_id"})[["uniprot_id"]]

    pe = pe.merge(selected, left_on="protein_id", right_on="uniprot_id", how="left", indicator=True)
    pe["group"] = pe["_merge"].map({"both": "known problematic", "left_only": "other"})

    labeled = pe.loc[pe["group"] == "known problematic", "p_vs_specific"]
    unlabeled = pe.loc[pe["group"] == "other", "p_vs_specific"]
    stat, pval = mannwhitneyu(labeled, unlabeled, alternative="two-sided")
    print(f"Mann-Whitney U: stat={stat:.1f} p={pval:.2e} (n_labeled={len(labeled)}, n_unlabeled={len(unlabeled)})")
    print(f"median p_vs_specific: known problematic={labeled.median():.3f} other={unlabeled.median():.3f}")

    fig, ax = plt.subplots(figsize=(4, 5))
    ax.boxplot([labeled, unlabeled], showfliers=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["known problematic", "other"])
    bracket(ax, 1, 2, pe["p_vs_specific"].quantile(0.98), pval)
    ax.set_ylabel("p_vs_specific")
    fig.tight_layout()
    fig.savefig(snakemake.output.enrichment, format="pdf")

    fig, ax = plt.subplots(figsize=(5, 4))
    for values, name in ((labeled, "known problematic"), (unlabeled, "other")):
        ax.hist(values, bins=40, range=(0, 1), density=True, alpha=0.5,
                label=f"{name} (n={len(values)})")
    ax.set_xlabel("p_vs_specific")
    ax.set_ylabel("density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(snakemake.output.histogram, format="pdf")


if __name__ == "__main__":
    main()
