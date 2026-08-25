import joblib
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import beta


def read_fasta_lengths(fasta_file):
    lengths = {}
    with open(fasta_file) as f:
        acc = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    lengths[acc] = length
                acc = line.split("|")[1]
                length = 0
            else:
                length += len(line.strip())
        if acc is not None:
            lengths[acc] = length
    return lengths


def interface_q(L, q, m, phi):
    mu = m.predict({"logL": np.log(L)})[0]
    a, b = mu * phi, (1 - mu) * phi
    return L * beta.ppf(q, a, b)


if __name__ == "__main__":
    interface_annotated_file = snakemake.input.interface_annotated
    fasta_file = snakemake.input.gene_fasta
    output_file = snakemake.output.model
    protein_lengths_file = snakemake.output.protein_lengths
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    df = pd.read_parquet(interface_annotated_file)
    print(f"Loaded {len(df)} annotated interactions", file=log, flush=True)

    lengths = read_fasta_lengths(fasta_file)
    df["length_bait"] = df["uniprot_id_bait"].map(lengths)
    df["length_prey"] = df["uniprot_id_prey"].map(lengths)
    df = df.dropna(subset=["length_bait", "length_prey"])

    print(f"{len(df)} interactions with known protein lengths", file=log, flush=True)

    df["min_length"] = df[["length_bait", "length_prey"]].min(axis=1)
    df["r"] = df["interface_residues"] / df["min_length"]

    df["logL"] = np.log(df.min_length)
    m = smf.glm("r ~ logL", df, family=sm.families.Binomial()).fit()
    phi = 1 / df.r.var()
    print(m.summary(), file=log, flush=True)

    fit_columns = [
        "uniprot_id_bait", "uniprot_id_prey",
        "length_bait", "length_prey", "min_length",
        "interface_residues", "r", "logL",
    ]
    df[fit_columns].to_csv(protein_lengths_file, index=False)

    # save model to file
    joblib.dump({"model": m, "phi": phi}, output_file)

    print("Done.", file=log, flush=True)
    log.close()
