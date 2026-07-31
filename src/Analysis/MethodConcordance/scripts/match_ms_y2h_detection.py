import pandas as pd


def load_detection(path, col_a, col_b, label):
    df = pd.read_parquet(path, columns=[col_a, col_b, "n_observed"])
    df = pd.DataFrame({
        "protein_a": df[[col_a, col_b]].min(axis=1),
        "protein_b": df[[col_a, col_b]].max(axis=1),
        f"{label}_detection": df["n_observed"] != 0,
    })
    return df.groupby(["protein_a", "protein_b"], as_index=False)[f"{label}_detection"].any()


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    col_a = snakemake.params.col_a
    col_b = snakemake.params.col_b

    ms_df = load_detection(snakemake.input.pod_ms, col_a, col_b, "ms")
    y2h_df = load_detection(snakemake.input.pod_y2h, col_a, col_b, "y2h")

    detection_df = ms_df.merge(y2h_df, on=["protein_a", "protein_b"], how="inner")
    detection_df.to_parquet(snakemake.output.detection, index=False)

    print(
        f"Matched pairs: {len(detection_df)} (MS pairs: {len(ms_df)}, Y2H pairs: {len(y2h_df)})",
        file=log, flush=True
    )
    log.close()
