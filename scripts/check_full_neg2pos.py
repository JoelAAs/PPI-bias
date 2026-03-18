
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

missing_baits_in_negative = positive_bait_prey_df[
    ~positive_bait_prey_df["bait"].isin(negative_bait_prey_df["bait"])
    ]
missing_prey_in_negative = positive_bait_prey_df[
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

    bait_target = {gene: round(degree * scale) for gene, degree in bait_target.items()}
    prey_target = {gene: round(degree * scale) for gene, degree in prey_target.items()}

    return bait_target, prey_target


def build_flow_graph(graph, target_bait, target_prey):
    F = nx.DiGraph()

    for node in graph.nodes():
        F.add_edge("source", ("bait", node), capacity=target_bait.get(node, 0).numerator)
        F.add_edge(("prey", node), "sink", capacity=target_prey.get(node, 0).numerator)

    for bait, prey in graph.edges():
        F.add_edge(("bait", bait), ("prey", prey), capacity=1)

    return F

for scaling in range(1, 2, 0.1):
    target_in, target_out = get_scaled_targets(positive_diG, scaling)
    print(f"Trying a subset where {pai} : 1")
    testF = build_flow_graph(negative_diG, target_in, target_out)
    flow_value, flow_dict = nx.maximum_flow(testF, "source", "sink")
    percent_output = round(flow_value / sum(target_in.values()).numerator * 100)
