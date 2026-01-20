import re

from ortools.sat.python import cp_model
import pandas as pd
import networkx as nx
import argparse

def read_edgelist(file):
    scaled_to = 1
    with open(file, "r") as f:
        for line in f:
            if line[0] == "#":
                match = re.search("Scaled: ([0-9]+) : 1", line)
                if match:
                    scaled_to = int(match.group(1))

    df = pd.read_csv(file, sep="\t", header=None, comment='#')
    df = df.rename({0: "bait", 1:"prey"}, axis=1)
    return df, scaled_to


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    parser.add_argument("--accepted_error", type=int, default=2, help="Path to output csv file")
    parser.add_argument("--threads", type=int, default=8, help="Path to output csv file")

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data
    balanced_positive = args.balanced_positive
    balanced_negative = args.balanced_negative
    missmatch_allowed = args.accepted_error
    threads = args.threads

    positive_df, _ = read_edgelist(positive_data)
    negative_df, scale = read_edgelist(negative_data)

    all_baits = set(positive_df["bait"]) | set(negative_df["bait"])
    all_preys = set(positive_df["prey"]) | set(negative_df["prey"])

    bait_int_idx = {b: i for i, b in enumerate(all_baits)}
    prey_int_idx = {p: i for i, p in enumerate(all_preys)}

    positive_diG = nx.from_pandas_edgelist(
        positive_df, "bait", "prey", create_using=nx.DiGraph())
    negative_diG = nx.from_pandas_edgelist(
        negative_df, "bait", "prey", create_using=nx.DiGraph())


    positive_edges = list(positive_diG.edges(data=True))
    negative_edges = list(negative_diG.edges(data=True))

    positive_edge_list_df = pd.DataFrame(positive_edges)
    positive_edge_list_df.columns = ["bait", "prey", "data"]

    negative_edge_list_df = pd.DataFrame(negative_edges)
    negative_edge_list_df.columns = ["bait", "prey", "data"]

    positive_edges_bait_ind = {}
    positive_edges_prey_ind = {}

    negative_edges_bait_ind = {}
    negative_edges_prey_ind = {}

    print("Indexing baits ...")
    for j, bait in enumerate(bait_int_idx):  # idx unused
        #print(f"{j}/{len(prey_int_idx)} ")
        negative_edges_bait_ind[bait] = negative_edge_list_df[negative_edge_list_df["bait"] == bait].index.to_list()
        positive_edges_bait_ind[bait] = positive_edge_list_df[positive_edge_list_df["bait"] == bait].index.to_list()

    print("Indexing prey ...")
    for j, prey in enumerate(prey_int_idx):
        print(f"{j}/{len(prey_int_idx)} ")
        negative_edges_prey_ind[prey]= negative_edge_list_df[negative_edge_list_df["prey"] == prey].index.to_list()
        positive_edges_prey_ind[prey]= positive_edge_list_df[positive_edge_list_df["prey"] == prey].index.to_list()

    print("Setting up model ...")
    model = cp_model.CpModel()
    xn = [model.NewBoolVar(f"xn_{i}") for i in range(len(negative_edges))]
    xp = [model.NewBoolVar(f"xp_{i}") for i in range(len(positive_edges))]

    K = model.NewIntVar(0, len(positive_edges) + len(negative_edges), "k")
    model.add(K == sum(xn) + sum(xp))
    bait_error = [model.NewIntVar(0, len(negative_edges), f"be_{i}")
                  for i in range(len(bait_int_idx))]
    prey_error = [model.NewIntVar(0, len(positive_edges), f"pe_{i}")
                  for i in range(len(prey_int_idx))]


    print("Setting up bait constraints ...")
    for bait, idx in bait_int_idx.items():
        s_n = sum(xn[i] for i in negative_edges_bait_ind[bait])
        s_p = sum(xp[i] for i in positive_edges_bait_ind[bait])
        diff = s_n - s_p*scale
        model.Add(bait_error[idx] >= diff - missmatch_allowed)
        model.Add(bait_error[idx] >= -diff - missmatch_allowed)
        model.Add(bait_error[idx] >= 0)
    # prey constraints

    print("Setting up prey constraints ...")
    for prey, idx in prey_int_idx.items():
        s_n = sum(xn[i] for i in negative_edges_prey_ind[prey])
        s_p = sum(xp[i] for i in positive_edges_prey_ind[prey])
        diff = s_n - s_p*scale
        model.Add(prey_error[idx] >= diff - missmatch_allowed)
        model.Add(prey_error[idx] >= -diff - missmatch_allowed)
        model.Add(prey_error[idx] >= 0)

    print("Minimizing degree delta and edges removed ...")
    model.Minimize(sum(bait_error) + sum(prey_error) - K*1.5)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = threads
    status = solver.Solve(model)
    print("Finished optimization")

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        total_mismatched_bait = 0
        total_mismatched_prey = 0
        n_negative_edges_selected = sum([solver.Value(xn[i]) for i, _ in enumerate(negative_edges)])
        n_positive_edges_selected = sum([solver.Value(xp[i]) for i, _ in enumerate(positive_edges)])

        for bait, idx in bait_int_idx.items():
            degree_negative = sum(solver.Value(xn[i]) for i, (u, v, d) in enumerate(negative_edges) if u == bait)
            degree_positive = sum(solver.Value(xp[i]) for i, (u, v, d) in enumerate(positive_edges) if u == bait)
            delta = abs(degree_negative - degree_positive*scale)
            total_mismatched_bait += delta

        for prey, idx in prey_int_idx.items():
            degree_negative = sum(solver.Value(xn[i]) for i, (u, v, d) in enumerate(negative_edges) if v == prey)
            degree_positive = sum(solver.Value(xp[i]) for i, (u, v, d) in enumerate(positive_edges) if v == prey)
            delta = abs(degree_negative - degree_positive*scale)
            total_mismatched_prey += delta

        print(f"Summed, bait-degree missmatch: {total_mismatched_bait}")
        print(f"Summed, Prey-degree missmatch: {total_mismatched_prey}")
        print(f"Negative edges retained: {n_negative_edges_selected}/{len(negative_edges)}")
        print(f"Positive edges retained: {n_positive_edges_selected}/{len(positive_edges)}")

        with open(balanced_negative, "w") as w:
            w.write(f"bait\tprey\n")
            for i, (u, v, d) in enumerate(negative_edges):
                if solver.Value(xn[i]) == 1:
                    w.write(f"{u}\t{v}\n")

        with open(balanced_positive, "w") as w:
            w.write(f"bait\tprey\n")
            for i, (u, v, d) in enumerate(positive_edges):
                if solver.Value(xp[i]) == 1:
                    w.write(f"{u}\t{v}\n")

    else:
        raise  ValueError("No feasible solution, Solver status:", status)
