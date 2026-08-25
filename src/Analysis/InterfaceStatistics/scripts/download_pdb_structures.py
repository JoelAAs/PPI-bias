import pandas as pd
from get_interfaces import download_pdb_structures, make_session


if __name__ == "__main__":
    matched_pairs_file = snakemake.input.matched_pairs
    pdb_folder = snakemake.output.pdb_folder
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    pp_df = pd.read_parquet(matched_pairs_file)
    print(f"Loaded {len(pp_df)} matched pairs", file=log, flush=True)

    with make_session() as session:
        download_pdb_structures(pp_df, session, pdb_folder, log)

    print("Done.", file=log, flush=True)
    log.close()
