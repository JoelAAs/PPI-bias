from ortools.sat.python import cp_model
import pandas as pd
import networkx as nx
import argparse

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

    positive_df = pd.read_csv(positive_data, sep="\t")
    negative_df = pd.read_csv(negative_data, sep="\t")

    positive_df = positive_df.iloc[:,0:2]
    negative_df = negative_df.iloc[:,0:2]

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


    positive_edges = list(positive_diG.edges(data=True))
    negative_edges = list(negative_diG.edges(data=True))

    positive_edge_list_df = pd.DataFrame(positive_edges)
    positive_edge_list_df.columns = ["bait", "prey", "data"]

    negative_edge_list_df = pd.DataFrame(negative_edges)
    negative_edge_list_df.columns = ["bait", "prey", "data"]

    # Convert to dict of lists
    print("Indexing edges ...", flush=True)
    neg_groups = negative_edge_list_df.groupby("bait").groups
    pos_groups = positive_edge_list_df.groupby("bait").groups

    negative_edges_bait_ind = {k: list(v) for k, v in neg_groups.items()}
    positive_edges_bait_ind = {k: list(v) for k, v in pos_groups.items()}
    
    print("Setting up model ...", flush=True)
    model = cp_model.CpModel()
    xn = [model.NewBoolVar(f"xn_{i}") for i in range(len(negative_edges))]
    xp = [model.NewBoolVar(f"xp_{i}") for i in range(len(positive_edges))]

    Kp = model.NewIntVar(0, len(positive_edges), "kp")
    Kn = model.NewIntVar(0, len(negative_edges), "kp")
    model.add(Kp == sum(xp))
    model.add(Kn == sum(xn))
    min_kept = len(positive_edges)
    model.add(Kp > int(len(positive_edges)*.9))
    model.add(Kn >= Kp)


    bait_error = [model.NewIntVar(0, len(negative_edges), f"be_{i}")
                  for i in range(len(bait_int_idx))]
    prey_error = [model.NewIntVar(0, len(positive_edges), f"pe_{i}")
                  for i in range(len(prey_int_idx))]


    print("Setting up bait constraints ...", flush=True)
    for bait, idx in bait_int_idx.items():
        s_n = sum(xn[i] for i in negative_edges_bait_ind[bait])
        s_p = sum(xp[i] for i in positive_edges_bait_ind[bait])
        diff = s_n - s_p
        model.add(bait_error[idx] >= diff - missmatch_allowed)
        model.add(bait_error[idx] >= -diff - missmatch_allowed)
        model.add(bait_error[idx] >= 0)
        # prey constraints

    print("Setting up prey constraints ...", flush=True)
    for prey, idx in prey_int_idx.items():
        s_n = sum(xn[i] for i in negative_edges_prey_ind[prey])
        s_p = sum(xp[i] for i in positive_edges_prey_ind[prey])
        diff = s_n - s_p
        model.add(prey_error[idx] >= diff - missmatch_allowed)
        model.add(prey_error[idx] >= -diff - missmatch_allowed)
        model.add(prey_error[idx] >= 0)

    print("Minimizing degree delta and edges removed ...", flush=True)
    model.Minimize(sum(bait_error) + sum(prey_error))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 600
    solver.parameters.log_search_progress = True
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
        print(f"Positive edges retained: {n_positive_edges_selected}/{len(positive_edges)}", flush=True)

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
