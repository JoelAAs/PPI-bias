import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pandas as pd


def endpoint_counts(table, filtered_table, id_col):
    is_pos = pc.greater_equal(table["n_observed"], 1).cast(pa.int64())
    pos_sub = table.select([id_col, "n_tested"]).append_column("is_pos", is_pos)
    pos_grouped = pos_sub.group_by(id_col).aggregate([
        ("n_tested", "sum"), ("is_pos", "sum"), (id_col, "count"),
    ]).rename_columns(["uniprot_id", "sum_tests", "deg_pos", "deg_tested"])

    is_neg = pc.equal(filtered_table["n_observed"], 0).cast(pa.int64())
    neg_sub = filtered_table.select([id_col]).append_column("is_neg", is_neg)
    neg_grouped = neg_sub.group_by(id_col).aggregate([("is_neg", "sum")]) \
        .rename_columns(["uniprot_id", "deg_neg"])

    return pos_grouped, neg_grouped


def bin_by_quantile(series, cutoff):
    hi = series.quantile(cutoff)
    lo = series[series != 0].quantile(1 - cutoff) # Remove all untested negative
    return np.where(series > hi, "high", np.where(series < lo, "low", "medium"))


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    min_tested = snakemake.params.hub_min_neg_tested
    hub_quantile_cutoff = snakemake.params.hub_quantile_cutoff
    col_a = snakemake.params.col_a
    col_b = snakemake.params.col_b

    dataset = ds.dataset(snakemake.input.pod)
    table = dataset.to_table(columns=[col_a, col_b, "n_tested", "n_observed"])
    print(f"Loaded {table.num_rows} pairs (undirected edge list) from {snakemake.input.pod}",
          file=log, flush=True)

    filtered_table = table.filter(pc.greater_equal(table["n_tested"], min_tested))
    print(f"{filtered_table.num_rows}/{table.num_rows} pairs have n_tested >= {min_tested} - "
          f"only these are eligible to count towards deg_neg", file=log, flush=True)

    a_pos, a_neg = endpoint_counts(table, filtered_table, col_a)
    b_pos, b_neg = endpoint_counts(table, filtered_table, col_b)

    pos_combined = pa.concat_tables([a_pos, b_pos]).group_by("uniprot_id").aggregate([
        ("sum_tests", "sum"), ("deg_pos", "sum"), ("deg_tested", "sum"),
    ]).rename_columns(["uniprot_id", "sum_tests", "deg_pos", "deg_tested"])

    neg_combined = pa.concat_tables([a_neg, b_neg]).group_by("uniprot_id").aggregate([
        ("deg_neg", "sum"),
    ]).rename_columns(["uniprot_id", "deg_neg"])

    degree_df = pos_combined.to_pandas().merge(neg_combined.to_pandas(), on="uniprot_id", how="left")
    degree_df["deg_neg"] = degree_df["deg_neg"].fillna(0).astype(int)
    print(f"{len(degree_df)} distinct proteins", file=log, flush=True)

    reference_counts = pd.read_csv(snakemake.input.reference_counts, sep="\t")
    degree_df = degree_df.merge(reference_counts, on="uniprot_id", how="left")
    n_unmatched = int(degree_df["n_pubmed"].isna().sum())
    degree_df["n_pubmed"] = degree_df["n_pubmed"].fillna(0).astype(int)
    print(f"Merged n_pubmed from {snakemake.input.reference_counts} "
          f"({n_unmatched} proteins with no UniProt reference_counts entry, set to 0)",
          file=log, flush=True)

    degree_df["deg_neg_bin"] = bin_by_quantile(degree_df["deg_neg"], hub_quantile_cutoff)
    degree_df["deg_pos_bin"] = bin_by_quantile(degree_df["deg_pos"], hub_quantile_cutoff)
    print(f"deg_neg_bin (cutoff={hub_quantile_cutoff}): "
          f"{degree_df['deg_neg_bin'].value_counts().to_dict()}", file=log, flush=True)
    print(f"deg_pos_bin (cutoff={hub_quantile_cutoff}): "
          f"{degree_df['deg_pos_bin'].value_counts().to_dict()}", file=log, flush=True)

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
