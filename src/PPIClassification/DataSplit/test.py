import numpy as np
import pandas as pd

positive_data = "work_folder/per_gene/subsets/validation/ms_limit_0.15_maxpos_pos.csv"
negative_data = "work_folder/per_gene/subsets/validation/ms_limit_3_poslim_0.15_maxpos_neg.csv"

positive_df = pd.read_csv(positive_data, sep="\t")
negative_df = pd.read_csv(negative_data, sep="\t")

positive_df = positive_df.iloc[:, 0:2]
negative_df = negative_df.iloc[:, 0:2]

positive_df.columns = ["bait", "prey"]
negative_df.columns = ["bait", "prey"]

positive_df = positive_df[
    (positive_df["bait"].isin(negative_df["bait"])) |
    (positive_df["prey"].isin(negative_df["prey"]))
    ]  # remove negative edges with no negative df representation

negative_df = negative_df[
    (negative_df["bait"].isin(positive_df["bait"])) |
    (negative_df["prey"].isin(positive_df["prey"]))
    ]

all_targets = {f"{g}_prey" for g in positive_df["prey"].unique()} | {f"{g}_bait" for g in positive_df["bait"].unique()}
gene_idx = {gene: index + 1 for index, gene in enumerate(all_targets)}

# gene_idx[0] is not in list
target_list = np.zeros(len(gene_idx) + 1,dtype=int)
for gene, c in positive_df.groupby("prey", as_index=False).size().values:
    target_list[gene_idx[f"{gene}_prey"]] = c

for gene, c in positive_df.groupby("bait", as_index=False).size().values:
    target_list[gene_idx[f"{gene}_bait"]] = c

target_list_save = target_list.copy()

negative_df["prey_idx"] = negative_df["prey"].apply(lambda gene: gene_idx.get(f"{gene}_prey", 0))
negative_df["bait_idx"] = negative_df["bait"].apply(lambda gene: gene_idx.get(f"{gene}_bait", 0))

chosen_edges_mask = np.zeros(negative_df.shape[0], dtype=int)
bait_idx = negative_df['bait_idx'].values
prey_idx = negative_df['prey_idx'].values

for i in range(negative_df.shape[0]):
    print(f"{i}/{negative_df.shape[0]}")
    bait, prey = bait_idx[i], prey_idx[i]
    if target_list[bait] > 0 and target_list[prey] > 0:
        chosen_edges_mask[i] = 1
        target_list[bait] -= 1
        target_list[prey] -= 1

def get_degree_change(value, idx):
    if idx == 0:
        return 0
    elif value == 0:
        return 1
    else:
        return -1

for i in range(negative_df.shape[0]):
    print(f"{i}/{negative_df.shape[0]}")
    if not chosen_edges_mask[i]:
        bait, prey = bait_idx[i], prey_idx[i]
        value_bait = target_list[bait]
        value_prey = target_list[prey]

        score_delta = get_degree_change(value_bait, bait) + get_degree_change(value_prey, prey)
        if score_delta < 0:
            chosen_edges_mask[i] = 1
            target_list[bait] -= 1
            target_list[prey] -= 1

print(f"{sum(target_list[1:])/sum(target_list_save)} % of edges missing")


import networkx as nx

G = nx.DiGraph()

# Add source/sink
G.add_node("s")
G.add_node("t")

# Add bait nodes
for bait in positive_df["bait"].unique():
    cap = target_list_save[gene_idx[f"{bait}_bait"]]
    G.add_edge("s", f"bait_{bait}", capacity=cap)

# Add prey nodes
for prey in positive_df["prey"].unique():
    cap = target_list_save[gene_idx[f"{prey}_prey"]]
    G.add_edge(f"prey_{prey}", "t", capacity=cap)

# Add negative edges
for _, row in negative_df.iterrows():
    G.add_edge(
        f"bait_{row['bait']}",
        f"prey_{row['prey']}",
        capacity=1
    )

flow_value, flow_dict = nx.maximum_flow(G, "s", "t")

bait_target = positive_df.groupby("bait").size().to_dict()
prey_target = positive_df.groupby("prey").size().to_dict()

constrained_baits = set(bait_target.keys())
constrained_preys = set(prey_target.keys())
for scale in range(10,0,-1):
    G = nx.DiGraph()

    SOURCE = "SOURCE"
    SINK = "SINK"

    G.add_node(SOURCE)
    G.add_node(SINK)

    INF = 10**12  # sufficiently large

    all_baits = set(negative_df["bait"].unique())

    for bait in all_baits:
        node_name = f"bait::{bait}"
        G.add_node(node_name)

        capacity = bait_target.get(bait, INF)
        G.add_edge(SOURCE, node_name, capacity=capacity*scale)

    all_preys = set(negative_df["prey"].unique())

    for prey in all_preys:
        node_name = f"prey::{prey}"
        G.add_node(node_name)

        capacity = prey_target.get(prey, INF)
        G.add_edge(node_name, SINK, capacity=capacity*scale)

    for _, row in negative_df.iterrows():
        bait_node = f"bait::{row['bait']}"
        prey_node = f"prey::{row['prey']}"
        G.add_edge(bait_node, prey_node, capacity=1)

    flow_value, flow_dict = nx.maximum_flow(G, SOURCE, SINK)
    print("Total flow:", 100*flow_value/(scale*positive_df.shape[0]), "% Scale:", scale )

chosen_edges = []
for _, row in negative_df.iterrows():
    bait_node = f"bait::{row['bait']}"
    prey_node = f"prey::{row['prey']}"

    if bait_node in flow_dict and prey_node in flow_dict[bait_node]:
        if flow_dict[bait_node][prey_node] > 0:
            chosen_edges.append((row["bait"], row["prey"]))


achieved_bait = {bait: 0 for bait in constrained_baits}
achieved_prey = {prey: 0 for prey in constrained_preys}

for bait, prey in chosen_edges:
    if bait in achieved_bait:
        achieved_bait[bait] += 1
    if prey in achieved_prey:
        achieved_prey[prey] += 1

print("Bait coverage:",
      sum(achieved_bait.values()) / sum(bait_target.values()))

print("Prey coverage:",
      sum(achieved_prey.values()) / sum(prey_target.values()))