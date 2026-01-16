import networkx as nx
import numpy as np
import pandas as pd
import argparse


def get_degree(diG):
    target_in_deg = diG.in_degree()
    target_out_deg = diG.out_degree()
    return target_in_deg, target_out_deg


def get_degree_list(nodes_ids, degree_view):
    degree_list = np.zeros(len(nodes_ids), dtype=int)
    for e, d in dict(degree_view).items():
        degree_list[e] = d
    return list(degree_list)


def get_degree_difference(targetG, sampledG):
    all_nodes = set(targetG.nodes()) | set(sampledG.nodes())
    degree_target = get_degree(targetG)
    degree_sampled = get_degree(sampledG)
    get_delta = lambda nodes, i, degree1, degree2: sum([
        np.abs(dict(degree1[i]).get(n, 0) - dict(degree2[i]).get(n, 0)) for n in nodes
    ])
    in_delta = get_delta(all_nodes, 0, degree_target, degree_sampled)
    out_delta = get_delta(all_nodes, 1, degree_target, degree_sampled)
    return in_delta + out_delta


def get_sampled_graph(in_degree, out_degree, sampling_G):
    s_diG = nx.DiGraph(nx.directed_configuration_model(
        in_degree, out_degree, sampling_G
    ))
    s_diG.remove_edges_from(nx.selfloop_edges(s_diG))
    return s_diG


def set_node_id(df, ids_dict):
    df.loc[:, "bait_id"] = df["bait"].map(ids_dict)
    df.loc[:, "prey_id"] = df["prey"].map(ids_dict)
    return df

def get_first_non_zero(degree_list, skip_single=True):
    for i, val in enumerate(degree_list):
        if val > int(skip_single):
            return i
    return None


def get_subsample_graph(targetG, sampledG, node_ids, size_setting="max", accepted_error=0.1):
    in_degree_target, out_degree_target = get_degree(targetG)
    in_degree_target_list = get_degree_list(node_ids.values(), in_degree_target)
    out_degree_target_list = get_degree_list(node_ids.values(), out_degree_target)
    n_pos_edges = len(targetG.edges())
    n_neg_edges = len(sampledG.edges())
    fraction_sampled_edges = n_pos_edges / n_neg_edges
    if size_setting == "max":
        for scaling_factor in range(100, 0, -1):
            scaled_in_degree_target_list = [
                round(d + d * fraction_sampled_edges * scaling_factor) for d in in_degree_target_list]
            scaled_out_degree_target_list = [
                round(d + d * fraction_sampled_edges * scaling_factor) for d in out_degree_target_list]
            degree_diff = sum(scaled_in_degree_target_list) - sum(scaled_out_degree_target_list)
            while degree_diff != 0:
                print(degree_diff)
                if degree_diff > 0:
                    scaled_in_degree_target_list[get_first_non_zero(scaled_in_degree_target_list)] -= 1
                else:
                    scaled_out_degree_target_list[get_first_non_zero(scaled_out_degree_target_list)] -= 1
            degree_diff = sum(scaled_in_degree_target_list) - sum(scaled_out_degree_target_list)

            s_diG = get_sampled_graph(scaled_in_degree_target_list, scaled_out_degree_target_list, sampledG)
            delta_degree = get_degree_difference(targetG, s_diG)
            if delta_degree < accepted_error * n_neg_edges:
                return s_diG, delta_degree, delta_degree / n_neg_edges
        return s_diG, delta_degree, delta_degree / n_neg_edges
    elif size_setting == "equal":
        s_diG = get_sampled_graph(in_degree_target_list, out_degree_target_list, sampledG)
        delta_degree = get_degree_difference(targetG, s_diG)
        return s_diG, delta_degree, delta_degree / n_neg_edges


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--balanced_negative", required=True, help="Path to output csv file")
    parser.add_argument("--balanced_positive", required=True, help="Path to output csv file")
    parser.add_argument("--size", default="max", help="Path to output csv file")
    parser.add_argument("--accepted_error", type=float, default=0.1, help="Path to output csv file")

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data
    balanced_positive = args.balanced_positive
    balanced_negative = args.balanced_negative
    size = args.size
    accepted_error = args.accepted_error

    positive_bait_prey_df = pd.read_csv(positive_data, sep="\t")
    negative_bait_prey_df = pd.read_csv(negative_data, sep="\t")

    positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
    negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

    positive_bait_prey_df.columns = ["bait", "prey"]
    negative_bait_prey_df.columns = ["bait", "prey"]

    baits = set(negative_bait_prey_df["bait"]) & set(positive_bait_prey_df["bait"])
    all_prey = set(negative_bait_prey_df["prey"]) & set(positive_bait_prey_df["prey"])

    negative_bp_df = negative_bait_prey_df[
        negative_bait_prey_df["bait"].isin(baits) & negative_bait_prey_df["prey"].isin(all_prey)].copy()
    positive_bp_df = positive_bait_prey_df[
        positive_bait_prey_df["bait"].isin(baits) & positive_bait_prey_df["prey"].isin(all_prey)].copy()

    node_ids = {gene_name: i for i, gene_name in enumerate(baits | all_prey)}
    id_to_gene = {i: gene for gene, i in node_ids.items()}
    positive_bp_df = set_node_id(positive_bp_df, node_ids)
    negative_bp_df = set_node_id(negative_bp_df, node_ids)
    positive_diG = nx.from_pandas_edgelist(
        positive_bp_df, "bait_id", "prey_id", create_using=nx.DiGraph()
    )
    negative_diG = nx.from_pandas_edgelist(
        negative_bp_df, "bait_id", "prey_id", create_using=nx.DiGraph()
    )

    sampledG, n_degree_dif, fraction_degree_diff = get_subsample_graph(
        positive_diG,
        negative_diG,
        node_ids,
        size_setting=size,
        accepted_error=accepted_error)

    print(f"degree differance {n_degree_dif}, at a percent of {fraction_degree_diff}")
    with open(balanced_negative, "w") as w:
        w.write("bait\tprey\n")
        for edge in sampledG.edges():
            bait, prey = edge
            w.write(f"{id_to_gene[bait]}\t{id_to_gene[prey]}\n")

    positive_bp_df[["bait","prey"]].to_csv(balanced_positive, sep="\t", index=False)