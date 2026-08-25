from get_interfaces import (
    load_protein_pair_file,
    build_entry_accessions,
    match_protein_pairs_to_pdbs,
    make_session,
)


if __name__ == "__main__":
    pod_file = snakemake.input.pod_file
    sifts_path = snakemake.input.sifts_file
    output_file = snakemake.output.matched_pairs
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    human_proteins, pp_df = load_protein_pair_file(pod_file)
    print(f"Loaded {len(pp_df)} pairs over {len(human_proteins)} proteins", file=log, flush=True)
    ent, ent_chains = build_entry_accessions(sifts_path, human_proteins, log)

    with make_session() as session:
        pp_df = match_protein_pairs_to_pdbs(pp_df, ent, session, log)

    pp_df.to_parquet(output_file, index=False)
    print("Done.", file=log, flush=True)
    log.close()
