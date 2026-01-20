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

    from ortools.sat.python import cp_model
    import pandas as pd
    import networkx as nx
    import numpy as np

    # ----------------------------
    # Load data
    # ----------------------------
    positive_data = "work_folder/per_gene/subsets/train/ms_pos.csv"
    negative_data = "work_folder/per_gene/subsets/train/ms_neg.csv"

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


    for bait, idx in bait_int_idx.items():
        s_n = sum(xn[i] for i,(u,v,d) in enumerate(negative_edges) if u == bait)
        s_p = sum(xp[i] for i,(u,v,d) in enumerate(positive_edges) if u == bait)
        diff = s_n - s_p
        model.Add(bait_error[idx] >= diff)
        model.Add(bait_error[idx] >= -diff)

    # prey constraints
    for prey, idx in prey_int_idx.items():
        s_n = sum(xn[i] for i, (u, v, d) in enumerate(negative_edges) if v == prey)
        s_p = sum(xp[i] for i, (u, v, d) in enumerate(positive_edges) if v == prey)
        diff =  s_n - s_p
        model.Add(prey_error[idx] >= diff)
        model.Add(prey_error[idx] >= -diff)


    model.minimize(sum(prey_error) + sum(bait_error) - K)
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    # ----------------------------
    # Inspect solution
    # ----------------------------

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        for bait, idx in bait_int_idx.items():
            degree_negative = sum(solver.Value(xn[i]) for i, (u, v, d) in enumerate(negative_edges) if u == bait)
            degree_positive = sum(solver.Value(xp[i]) for i, (u, v, d) in enumerate(positive_edges) if u == bait)
            delta = abs(degree_negative - degree_positive)
            print(bait,"bp:", degree_positive,"bn:", degree_negative, "delta:", delta)
    else:
        print("No feasible solution found! Solver status:", status)

