import numpy as np
import pandas as pd


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    
    l_q = snakemake.params.low_quantile

    lof = pd.read_csv(snakemake.input.gnomad_lof_metrics, sep="\t", compression="gzip")
    lof = lof[["gene", "oe_lof_upper"]].dropna().drop_duplicates(subset="gene")
    print(f"{len(lof)} genes with an oe_lof_upper value", file=log, flush=True)

    threshold = lof["oe_lof_upper"].quantile(l_q)
    lof["category"] = np.where(lof["oe_lof_upper"] <= threshold, "low", "high")
    print(
        f"oe_lof_upper {l_q:.0%} threshold: {threshold:.4f} "
        f"({(lof['category'] == 'low').sum()} low, {(lof['category'] == 'high').sum()} high)",
        file=log, flush=True,
    )

    symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t")
    category_df = symbol_map.merge(lof, left_on="hgnc_symbol", right_on="gene", how="inner")
    print(
        f"{category_df['uniprot_id'].nunique()} uniprot IDs mapped to an oe_lof_upper category",
        file=log, flush=True,
    )

    category_df[["uniprot_id", "gene", "oe_lof_upper", "category"]].to_csv(
        snakemake.output.oe_lof_category, sep="\t", index=False
    )
    log.close()
