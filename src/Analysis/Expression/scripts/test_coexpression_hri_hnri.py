import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def attach_coexpression(df, uniprot_to_row, corr):
    row_a = df["bait"].map(uniprot_to_row)
    row_b = df["prey"].map(uniprot_to_row)
    valid = row_a.notna() & row_b.notna()
    n_dropped = (~valid).sum()
    out = df.loc[valid].copy()
    out["coexpr"] = corr[row_a[valid].to_numpy(dtype=int), row_b[valid].to_numpy(dtype=int)]
    return out, n_dropped


def bin_or_between_groups(combined, bin_value):
    in_bin = combined["expression_bin"] == bin_value
    pos = combined["pos"] == 1
    table = [
        [(in_bin & pos).sum(), (~in_bin & pos).sum()],
        [(in_bin & ~pos).sum(), (~in_bin & ~pos).sum()],
    ]
    odds_ratio, p_value = fisher_exact(table)
    return odds_ratio, p_value


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    gene_index = pd.read_csv(snakemake.input.gene_index, sep="\t")
    uniprot_to_row = dict(zip(gene_index["uniprot_id"], gene_index["row_index"]))
    corr = np.load(snakemake.input.corr_npy)

    hrni = pd.read_csv(snakemake.input.balanced_negative, sep="\t", dtype={"bait": "string", "prey": "string"})
    hri = pd.read_csv(snakemake.input.balanced_positive, sep="\t", dtype={"bait": "string", "prey": "string"})
    hri, n_drop_hri = attach_coexpression(hri, uniprot_to_row, corr)
    hrni, n_drop_hrni = attach_coexpression(hrni, uniprot_to_row, corr)

    hri = hri.reset_index(drop=True)
    hrni = hrni.reset_index(drop=True)
    combined = pd.concat([hri.assign(pos=1), hrni.assign(pos=0)], ignore_index=True)

    corr_values = corr[~np.isnan(corr)]
    low_q = np.quantile(corr_values, 0.15)
    high_q = np.quantile(corr_values, 0.85)
    combined["expression_bin"] = np.where(
        combined["coexpr"] <= low_q, "low",
        np.where(combined["coexpr"] >= high_q, "high", "medium"),
    )

    or_high, p_high = bin_or_between_groups(combined, "high")
    or_low, p_low = bin_or_between_groups(combined, "low")
    log.write(f"OR high vs rest (pos vs neg): {or_high} (p={p_high})\n")
    log.write(f"OR low vs rest (pos vs neg): {or_low} (p={p_low})\n")

    with open(snakemake.output.expression_OR, "w") as w:
        msg = (
            f"OR\tp-value\tclass\n"
            f"{or_high}\t{p_high}\thigh\n"
            f"{or_low}\t{p_low}\tlow\n"
        )
        w.write(msg)


    combined.to_csv(snakemake.output.hrni_vs_hri, sep="\t", index=None)

