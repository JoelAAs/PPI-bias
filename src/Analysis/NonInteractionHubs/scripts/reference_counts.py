import pandas as pd

if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    df = pd.read_csv(snakemake.input.uniprot_annotation, sep="\t", dtype={"lit_pubmed_id": str})
    print(f"Loaded {len(df)} human proteins from {snakemake.input.uniprot_annotation}",
          file=log, flush=True)

    def count_pubmeds(cell):
        if pd.isna(cell):
            return 0
        return len({pid.strip() for pid in cell.split(";") if pid.strip()})

    counts = pd.DataFrame({
        "uniprot_id": df["uniprot_id"],
        "n_pubmed": df["lit_pubmed_id"].apply(count_pubmeds),
    })
    n_zero = int((counts["n_pubmed"] == 0).sum())
    print(f"{len(counts)} proteins, {n_zero} with no lit_pubmed_id entries in UniProt",
          file=log, flush=True)

    counts.to_csv(snakemake.output.reference_counts, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
