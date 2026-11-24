import numpy as np
import pandas as pd
from balance_naive_sets import drop_shared_pairs


def main():
    df_neg = pd.read_csv(snakemake.input.neg, sep="\t")
    df_pos = pd.read_csv(snakemake.input.pos, sep="\t")
    if snakemake.wildcards.pos_set == "hippie":
        df_pos = df_pos[df_pos["bait"] != df_pos["prey"]]
        df_pos = df_pos[~(df_pos["bait"].str.contains("-") | df_pos["prey"].str.contains("-"))]
        swap = (df_pos["bait"] > df_pos["prey"]).to_numpy()
        df_pos = pd.DataFrame({
            "bait": np.where(swap, df_pos["prey"], df_pos["bait"]),
            "prey": np.where(swap, df_pos["bait"], df_pos["prey"])
        })
        df_pos = df_pos.drop_duplicates()

    df_pos, df_neg = drop_shared_pairs(df_pos, df_neg)

    df_pos[["bait", "prey"]].rename(columns={
        "bait": "protein1", "prey": "protein2"
    }).to_csv(snakemake.output.pos, sep=",", index=False)
    df_neg[["bait", "prey"]].rename(columns={
        "bait": "protein1", "prey": "protein2"
    }).to_csv(snakemake.output.neg, sep=",", index=False)


if __name__ == "__main__":
    main()
