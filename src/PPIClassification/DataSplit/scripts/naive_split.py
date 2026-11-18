import sys
import numpy as np
import pandas as pd


def split(df, test_fraction, rng):
    proteins = np.array(sorted(set(df["bait"]) | set(df["prey"])))  # sorted: set order changes between runs
    test = rng.choice(proteins, size=int(len(proteins) * test_fraction), replace=False)
    train_df = df[~(df["bait"].isin(test)) & ~(df["prey"].isin(test))]
    test_df = df[(df["bait"].isin(test)) & (df["prey"].isin(test))]
    return train_df, test_df


def write_dataset(pos_df, neg_df, output_file):
    pd.concat([pos_df, neg_df], ignore_index=True).to_csv(output_file, sep="\t", index=False)


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    rng = np.random.default_rng(snakemake.params.seed)
    test_fraction = snakemake.params.test_fraction
    df = pd.read_csv(snakemake.input.pos, sep="\t", dtype={"bait": str, "prey": str})
    train_pos, test_pos = split(df, test_fraction, rng)
    train_pos.to_csv(snakemake.output.train, sep="\t", index=False)
    test_pos.to_csv(snakemake.output.test, sep="\t", index=False)

if __name__ == "__main__":
    main()
