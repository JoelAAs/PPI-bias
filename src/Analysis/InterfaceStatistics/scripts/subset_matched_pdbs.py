import pandas as pd


if __name__ == "__main__":
    all_matched_file = snakemake.input.all_matched_pairs
    pod_file = snakemake.input.pod_file
    output_file = snakemake.output.matched_pairs
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    all_matched = pd.read_parquet(all_matched_file)
    pod_df = pd.read_parquet(pod_file)

    pdb_lookup = all_matched.set_index(["uniprot_id_bait", "uniprot_id_prey"])["matching_pdb"]
    keys = pd.MultiIndex.from_arrays([pod_df["uniprot_id_bait"], pod_df["uniprot_id_prey"]])
    pod_df["matching_pdb"] = pdb_lookup.reindex(keys).to_numpy()
    pod_df = pod_df[pod_df["matching_pdb"].notna()].copy()

    print(f"{len(pod_df)}/{len(all_matched)} flat-matched pairs present in this dataset", file=log, flush=True)
    pod_df.to_parquet(output_file, index=False)
    print("Done.", file=log, flush=True)
    log.close()
