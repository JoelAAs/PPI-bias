
import networkx as nx
from fractions import Fraction
import pandas as pd
import argparse


negative_data = "work_folder/per_gene/subsets/ms_directional_full_1_neg.pq"
positive_data = "work_folder/per_gene/subsets/ms_directional_full_0.02_pos.pq"

positive_bait_prey_df = pd.read_parquet(positive_data)
negative_bait_prey_df = pd.read_parquet(negative_data)

positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

positive_bait_prey_df.columns = ["bait", "prey"]
negative_bait_prey_df.columns = ["bait", "prey"]

shared_baits = set(positive_bait_prey_df["bait"]) & set(negative_bait_prey_df["bait"])
shared_prey = set(positive_bait_prey_df["prey"]) & set(negative_bait_prey_df["prey"])

missing_baits_in_negative = negative_bait_prey_df[
    ~positive_bait_prey_df["bait"].isin(negative_bait_prey_df["bait"])
    ]
missing_prey_in_negative = negative_bait_prey_df[
    ~positive_bait_prey_df["prey"].isin(negative_bait_prey_df["prey"])
    ]

print("--------------------------EDGES DISCARDED---------------------------------")
print(f"{len(set(positive_bait_prey_df["bait"]) - shared_baits)} / {len(shared_baits)} of "
      f"baits are missing in the negative data")
print(f"This leads to {(~positive_bait_prey_df['bait'].isin(shared_baits)).sum()} "
      f"edges being discarded.")

print(f"{len(set(positive_bait_prey_df["prey"]) - shared_prey)} / {len(shared_prey)} of "
      f"baits are missing in the negative data")
print(f"This leads to {(~positive_bait_prey_df['prey'].isin(shared_prey)).sum()} "
      f"edges being discarded.")


positive_bait_prey_df[positive_bait_prey_df["bait"].isin(shared_baits) & positive_bait_prey_df["prey"].isin(shared_prey)]
negative_bait_prey_df[negative_bait_prey_df["bait"].isin(shared_baits) & negative_bait_prey_df["prey"].isin(shared_prey)]

positive_diG = nx.from_pandas_edgelist(positive_bait_prey_df, "bait", "prey", create_using=nx.DiGraph())
negative_diG = nx.from_pandas_edgelist(negative_bait_prey_df, "bait", "prey", create_using=nx.DiGraph())
success = False

def get_scaled_targets(graph, scale):
    bait_target = dict(graph.out_degree())
    prey_target = dict(graph.in_degree())

    bait_target = {gene: round(degree * scaling) for gene, degree in bait_target.items()}
    prey_target = {gene: round(degree * scaling) for gene, degree in prey_target.items()}

    return bait_target, prey_target

for scaling in range(1, 2, 0.1):
    target_in, target_out, success = get_scaled_targets(positive_diG, scaling)
    if success:
        print(f"Trying a subset where {pai} : 1")
        testF = build_flow_graph(negative_diG, target_in, target_out)
        flow_value, flow_dict = nx.maximum_flow(testF, "source", "sink")
        percent_output = round(flow_value / sum(target_in.values()).numerator * 100)

        print(f"Flow value: {flow_value}, that being {percent_output} % of scaled degree")

        min_target_ppis = sum(target_in.values()) / pai
        min_ppi_target = min_target_ppis*.8 < flow_value < min_target_ppis*1.2

        save=False
        if args.subset == "test":
            if min_ppi_target or pai == 1:
                save = True
        elif percent_output > min_max_flow or pai == 1 or min_ppi_target: # Fix this one later
            save = True

        if save:
            S = nx.DiGraph()
            S.add_nodes_from(negative_diG.nodes())

            for u, v in negative_diG.edges():
                if flow_dict.get(("out", u), {}).get(("in", v), 0) == 1:
                    S.add_edge(u, v)

            with open(max_flow_negative, "w") as w:
                w.write(f"#Scaled: {pai.numerator} : 1\n")
                for u, v in S.edges():
                    w.write(f"{u}\t{v}\n")

            with open(max_flow_positive, "w") as w:
                for u, v in positive_diG.edges():
                    w.write(f"{u}\t{v}\n")
            success = True
            break
if not success:
    raise ValueError(f"No possible subset with flow > {min_max_flow} %")

