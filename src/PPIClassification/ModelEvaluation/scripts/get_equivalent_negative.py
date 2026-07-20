import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "DataSplit" / "scripts"))
from sample_balance_multi_network_functions import degree_balace_edges


def main():
    sel_neg_df = pd.read_csv(snakemake.input.selected_negative, sep="\t", dtype={"bait": "string", "prey": "string"})
    pos_df = pd.read_csv(snakemake.input.test_pos, sep="\t", dtype={"bait": "string", "prey": "string"})
    network_type = snakemake.wildcards.network_type
    if network_type == "directional":
        directed = True
    elif network_type == "undirectional":
        directed = False
    else:
        raise ValueError(f"{network_type} is not a valid network type.")
    
    all_nodes = list(set(sel_neg_df["bait"]) | set(sel_neg_df["prey"]) )
    pos_set = set(zip(pos_df["bait"], pos_df["prey"])) | set(zip(sel_neg_df["bait"], sel_neg_df["prey"]))
    comp_df = pd.DataFrame(
        [(u, v) for u in all_nodes for v in all_nodes if u != v and (u, v) not in pos_set],
        columns=["bait", "prey"],
    )
    if not directed:
        comp_df["id"] = comp_df.apply(lambda row: "-".join(sorted(row)), axis=1)
        comp_df.drop_duplicates(inplace=True)
        del comp_df["id"] # no duplicates for faster loading later

    _, selected_neg, _ = degree_balace_edges(sel_neg_df, comp_df, min_flow=0.80, directed=directed)

    selected_neg.to_csv(snakemake.output.non_obs, sep="\t", index=False)


if __name__ == "__main__":
    main()
    