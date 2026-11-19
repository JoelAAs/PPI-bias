import sys
import numpy as np
import pandas as pd
from random_negative import get_no_edges



def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    directed = {"directional": True, "undirectional": False}[snakemake.wildcards.network_type]

    df = pd.read_csv(snakemake.input.balanced, sep="\t")
    df_pos = df[df["interaction"]]
    df_neg = df[~df["interaction"]]
    selected_negative_edges = set(zip(df_neg["bait"], df_neg["prey"]))
        
    no_neg_df = get_no_edges(df_pos, directed=directed, edges_to_avioid=selected_negative_edges)
    print(f"Observed negative: {df_neg.shape[0]}, selected NO negative: {no_neg_df.shape[0]}", flush=True)

    balanced_no = pd.concat(
        [df_pos, no_neg_df.assign(interaction=False)],
        ignore_index=True
    )
    balanced_no[["bait", "prey", "interaction"]].to_csv(snakemake.output.balanced_no, sep="\t", index=False)


if __name__ == "__main__":
    main()
