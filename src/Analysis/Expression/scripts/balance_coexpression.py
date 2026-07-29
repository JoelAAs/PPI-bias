import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "PPIClassification" / "DataSplit" / "scripts"))
from sample_balance_multi_network_functions import degree_balace_edges


def main():
    pos_df = pd.read_csv(
        snakemake.input.max_positive, sep="\t", header=0,
        names=["bait", "prey"], dtype={"bait": "string", "prey": "string"}
    )
    neg_df = pd.read_csv(
        snakemake.input.max_negative, sep="\t", header=0,
        names=["bait", "prey"], dtype={"bait": "string", "prey": "string"}
    )
    
    expression_uniprotids = set(pd.read_csv(snakemake.input.gene_index, sep="\t")["uniprot_id"])
    neg_df = neg_df[
        (neg_df["bait"].isin(expression_uniprotids)) &
        (neg_df["prey"].isin(expression_uniprotids))
    ]
    pos_df = pos_df[
        (pos_df["bait"].isin(expression_uniprotids)) &
        (pos_df["prey"].isin(expression_uniprotids))
    ]

    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")

    selected_pos, selected_neg, _ = degree_balace_edges(pos_df, neg_df, min_flow=0.95, directed=directed)

    selected_pos.to_csv(snakemake.output.balanced_positive, sep="\t", index=False)
    selected_neg.to_csv(snakemake.output.balanced_negative, sep="\t", index=False)


if __name__ == "__main__":
    main()
    