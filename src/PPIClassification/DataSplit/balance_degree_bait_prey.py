import math
import random
from fractions import Fraction

import numpy as np
import networkx as nx
import pandas as pd
from ortools.sat.python import cp_model


def get_degrees(bait_idx, prey_idx, G):
    target_baits = np.zeros(len(bait_idx), dtype=int)
    target_preys = np.zeros(len(prey_idx), dtype=int)
    for bait, prey in G.edges():
        target_baits[bait_idx[bait]] += 1
        target_preys[prey_idx[prey]] += 1

    return target_baits, target_preys

def get_error(target_degrees, s_degree):
    return np.abs(target_degrees[0] - s_degree[0]).sum() + np.abs(target_degrees[0] - s_degree[0]).sum()

def scaling_fractions(min_edges, max_edges):
    return [
        Fraction(i, min_edges) for i in range(max_edges, 1, -1)
    ]


def get_possible_scaling_factors(targetG, scaling):
    target_degree_in, target_degree_out = get_degrees(targetG)
    target_in_scaled = {}
    target_out_scaled = {}
    feasible = True
    for node in targetG.nodes():
        target_in_scaled[node] = target_degree_in[node] * scaling
        target_out_scaled[node] = target_degree_out[node] * scaling
        if target_in_scaled[node].denominator != 1 or target_out_scaled[node].denominator != 1:
            feasible = False
            break
    return target_in_scaled, target_out_scaled, feasible

def update_drop(s_degree, drop_idx):


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--selected_ppis", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    parser.add_argument("--subtractive", default=False, action="store_true", help="Path to output csv file")
    parser.add_argument("--size", default="max", help="Path to output csv file")
    parser.add_argument("--accepted_error", type=float, default=0.1, help="Path to output csv file")

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data

    selected_ppi_file = args.selected_ppis

    balanced_positive = args.balanced_positive
    balanced_negative = args.balanced_negative

    subtractive = args.subtractive
    size = args.size
    accepted_error = args.accepted_error





    positive_data = "work_folder/per_gene/subsets/train/ms_pos.csv"
    negative_data = "work_folder/per_gene/subsets/train/ms_neg.csv"
    positive_bait_prey_df = pd.read_csv(positive_data, sep="\t")
    negative_bait_prey_df = pd.read_csv(negative_data, sep="\t")

    positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
    negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

    positive_bait_prey_df.columns = ["bait", "prey"]
    negative_bait_prey_df.columns = ["bait", "prey"]

    all_baits = set(positive_bait_prey_df["bait"]) | set(negative_bait_prey_df["bait"])
    all_prey = set(positive_bait_prey_df["prey"]) | set(negative_bait_prey_df["prey"])

    bait_int_idx = {bait: i for i, bait in enumerate(all_baits)}
    prey_int_idx = {prey: i for i, prey in enumerate(all_prey)}

    positive_diG = nx.from_pandas_edgelist(
        positive_bait_prey_df, "bait", "prey", create_using=nx.DiGraph())
    negative_diG = nx.from_pandas_edgelist(
        negative_bait_prey_df, "bait", "prey", create_using=nx.DiGraph())

    negative_edges = list(negative_diG.edges(data=True))
    target_degrees = get_degrees(bait_int_idx, prey_int_idx, positive_diG)


    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x_{i}") for i in range(len(negative_edges))]
    bait_error = [model.NewIntVar(0, len(negative_edges), f"be_{i}")
                  for i in range(len(bait_int_idx))]
    prey_error = [model.NewIntVar(0, len(negative_edges), f"pe_{i}")
                  for i in range(len(prey_int_idx))]


    pa = scaling_fractions(len(positive_diG.edges()), len(negative_diG.edges()))
    for pai in pa:
        target_in, target_out, success = get_possible_scaling_factors(positive_diG, pai)


    e = 0.1
    upper_bound_bait = list(map(target_degrees[0], lambda x: math.ceil(x*(1+e))))
    lower_bound_bait = list(map(target_degrees[0], lambda x: math.ceil(x*(1+e))))

    for bait, idx in bait_int_idx.items():
        s = sum(x[i] for i,(u,v,d) in enumerate(negative_edges) if u == bait)
        model.Add(bait_error[idx] >= s - target_degrees[0][idx])
        model.Add(bait_error[idx] >= target_degrees[0][idx] - s)

    # prey constraints
    for prey, idx in prey_int_idx.items():
        s = sum(x[i] for i,(u,v,d) in enumerate(negative_edges) if v == prey)
        model.Add(prey_error[idx] >= s - target_degrees[1][idx])
        model.Add(prey_error[idx] >= target_degrees[1][idx] - s)

    #model.Minimize(sum(bait_error) + sum(prey_error))
    model.maximize()
    solver = cp_model.CpSolver()
    solver.Solve(model)
    selected_indices = [i for i in range(len(x)) if solver.Value(x[i]) == 1]
    subset_edges = [negative_edges[i] for i in selected_indices]


