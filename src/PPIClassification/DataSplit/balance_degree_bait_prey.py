from ortools.sat.python import cp_model
import pandas as pd
import networkx as nx
import numpy as np

def get_degrees(bait_idx, prey_idx, G):
    target_baits = np.zeros(len(bait_idx), dtype=int)
    target_preys = np.zeros(len(prey_idx), dtype=int)
    for bait, prey in G.edges():
        target_baits[bait_idx[bait]] += 1
        target_preys[prey_idx[prey]] += 1

    return target_baits, target_preys



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--selected_ppis", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    parser.add_argument("--subtractive", default=False, action="store_true", help="Path to output csv file")
    parser.add_argument("--size", default="max", help="Path to output csv file")
    parser.add_argument("--accepted_error", type=int, default=2, help="Path to output csv file")

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data

    selected_ppi_file = args.selected_ppis

    balanced_positive = args.balanced_positive
    balanced_negative = args.balanced_negative

    accepted_error = args.accepted_error


    # positive_data = "work_folder/per_gene/subsets/train/ms_pos.csv"
    # negative_data = "work_folder/per_gene/subsets/train/ms_neg.csv"

    positive_df = pd.read_csv(positive_data, sep="\t")[["gene_name_bait", "gene_name_prey"]]
    negative_df = pd.read_csv(negative_data, sep="\t")[["gene_name_bait", "gene_name_prey"]]

    positive_df.columns = ["bait", "prey"]
    negative_df.columns = ["bait", "prey"]

    all_baits = set(positive_df["bait"]) | set(negative_df["bait"])
    all_preys = set(positive_df["prey"]) | set(negative_df["prey"])

    bait_int_idx = {b: i for i, b in enumerate(all_baits)}
    prey_int_idx = {p: i for i, p in enumerate(all_preys)}

    positive_diG = nx.from_pandas_edgelist(
        positive_df, "bait", "prey", create_using=nx.DiGraph())
    negative_diG = nx.from_pandas_edgelist(
        negative_df, "bait", "prey", create_using=nx.DiGraph())

    negative_edges = list(negative_diG.edges(data=True))
    positive_edges = list(positive_diG.edges(data=True))

    model = cp_model.CpModel()
    xn = [model.NewBoolVar(f"xn_{i}") for i in range(len(negative_edges))]
    xp = [model.NewBoolVar(f"xp_{i}") for i in range(len(positive_edges))]

    K = model.NewIntVar(0, len(positive_edges) + len(negative_edges), "k")
    model.add(K == sum(xn) + sum(xp))
    bait_error = [model.NewIntVar(0, len(negative_edges), f"be_{i}")
                  for i in range(len(bait_int_idx))]
    prey_error = [model.NewIntVar(0, len(positive_edges), f"pe_{i}")
                  for i in range(len(prey_int_idx))]


    missmatch_allowed = accepted_error
    for bait, idx in bait_int_idx.items():
        s_n = sum(xn[i] for i,(u,v,d) in enumerate(negative_edges) if u == bait)
        s_p = sum(xp[i] for i,(u,v,d) in enumerate(positive_edges) if u == bait)
        diff = s_n - s_p
        model.Add(bait_error[idx] >= diff - missmatch_allowed)
        model.Add(bait_error[idx] >= -diff - missmatch_allowed)
        model.Add(prey_error[idx] >= 0)
    # prey constraints
    for prey, idx in prey_int_idx.items():
        s_n = sum(xn[i] for i, (u, v, d) in enumerate(negative_edges) if v == prey)
        s_p = sum(xp[i] for i, (u, v, d) in enumerate(positive_edges) if v == prey)
        diff =  s_n - s_p
        model.Add(prey_error[idx] >= diff - missmatch_allowed)
        model.Add(prey_error[idx] >= -diff - missmatch_allowed)
        model.Add(prey_error[idx] >= 0)


    model.minimize(sum(bait_error) + sum(prey_error) - K)
    solver = cp_model.CpSolver()
    status = solver.Solve(model)


    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        total_mismatched_bait = 0
        total_mismatched_prey = 0
        n_negative_edges_selected = sum([solver.Value(xn[i]) for i, _ in enumerate(negative_edges)])
        n_positive_edges_selected = sum([solver.Value(xp[i]) for i, _ in enumerate(positive_edges)])

        for bait, idx in bait_int_idx.items():
            degree_negative = sum(solver.Value(xn[i]) for i, (u, v, d) in enumerate(negative_edges) if u == bait)
            degree_positive = sum(solver.Value(xp[i]) for i, (u, v, d) in enumerate(positive_edges) if u == bait)
            delta = abs(degree_negative - degree_positive)
            total_mismatched_bait += delta
        for prey, idx in prey_int_idx.items():
            degree_negative = sum(solver.Value(xn[i]) for i, (u, v, d) in enumerate(negative_edges) if v == prey)
            degree_positive = sum(solver.Value(xp[i]) for i, (u, v, d) in enumerate(positive_edges) if v == prey)
            delta = abs(degree_negative - degree_positive)
            total_mismatched_prey += delta

        print(f"Summed, bait-degree missmatch: {total_mismatched_bait}")
        print(f"Summed, Prey-degree missmatch: {total_mismatched_prey}")
        print(f"Negative edges retained: {n_negative_edges_selected}/{len(negative_edges)}")
        print(f"Positive edges retained: {n_positive_edges_selected}/{len(positive_edges)}")

        with open("balanced_negative_edgelist_test.csv", "w") as w:
            w.write(f"bait\tprey\n")
            for i, (u, v, d) in enumerate(negative_edges):
                if solver.Value(xn[i]) == 1:
                    w.write(f"{u}\t{v}\n")

        with open("balanced_positive_edgelist_test.csv", "w") as w:
            w.write(f"bait\tprey\n")
            for i, (u, v, d) in enumerate(positive_edges):
                if solver.Value(xp[i]) == 1:
                    w.write(f"{u}\t{v}\n")

    else:
        print("No feasible solution found! Solver status:", status)

