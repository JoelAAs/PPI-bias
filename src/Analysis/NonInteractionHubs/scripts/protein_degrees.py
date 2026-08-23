import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pandas as pd


def endpoint_counts(table, id_col):
    is_pos = pc.greater_equal(table["n_observed"], 1).cast(pa.int64())
    sub = table.select([id_col, "n_tested"]).append_column("is_pos", is_pos)

    degrees = sub.group_by(id_col).aggregate([
        ("n_tested", "sum"), ("is_pos", "sum"), (id_col, "count"),
    ]).rename_columns(["uniprot_id", "sum_tests", "deg_pos", "deg_tested"])
    degrees = degrees.append_column("deg_neg", pc.subtract(degrees["deg_tested"], degrees["deg_pos"]))
    
    return degrees


def bin_by_quantile(series, cutoff):
    hi = series.quantile(cutoff)
    lo = series[series != 0].quantile(1 - cutoff) # Remove all untested negative
    return np.where(series > hi, "high", np.where(series < lo, "low", "medium"))


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    
    col_a = snakemake.params.col_a
    col_b = snakemake.params.col_b
    int_limits = snakemake.params.interaction_hubs
    non_hub_q = snakemake.params.non_interaction_hub_q

    dataset = ds.dataset(snakemake.input.pod)
    table = dataset.to_table(columns=[col_a, col_b, "n_tested", "n_observed"])
    
    degree_a = endpoint_counts(table, col_a)
    degree_b = endpoint_counts(table, col_b)

    combined = pa.concat_tables([degree_a, degree_b]).group_by("uniprot_id").aggregate([
        ("sum_tests", "sum"), ("deg_pos", "sum"), ("deg_neg", "sum"), ("deg_tested", "sum")
    ]).rename_columns(["uniprot_id", "sum_tests", "deg_pos", "deg_neg", "deg_tested"])

    degree_df = combined.to_pandas()
    
    reference_counts = pd.read_csv(snakemake.input.reference_counts, sep="\t")
    degree_df = degree_df.merge(reference_counts, on="uniprot_id", how="left")
    degree_df["n_pubmed"] = degree_df["n_pubmed"].fillna(0).astype(int)

    degree_df["deg_pos_bin"] = np.where(
        degree_df["deg_pos"] <= int_limits[0], "low",
        np.where(degree_df["deg_pos"] > int_limits[1], "high", "medium")
    ) # pos degree on biology

    degree_df["deg_neg_bin"] = bin_by_quantile(degree_df["deg_neg"], non_hub_q) # degree based

    degree_df.to_parquet(snakemake.output.degree, index=False)

    bin_order = ["low", "medium", "high"]
    bin_summary = pd.DataFrame({
        "max_deg_neg": degree_df.groupby("deg_neg_bin")["deg_neg"].max(),
        "max_deg_pos": degree_df.groupby("deg_pos_bin")["deg_pos"].max(),
    }).reindex(bin_order)
    bin_summary.index.name = "bin"
    bin_summary.to_csv(snakemake.output.bin_summary, sep="\t")
    print(f"Bin summary (highest deg_neg/deg_pos per bin):\n{bin_summary}", file=log, flush=True)

    print("Done.", file=log, flush=True)
    log.close()
