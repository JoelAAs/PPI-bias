import networkx as nx
from fractions import Fraction
import pandas as pd
import argparse

def get_degrees(G):
    degree_in = dict(G.in_degree())
    degree_out = dict(G.out_degree())
    return degree_in, degree_out


def scaling_fractions(min_edges, max_edges):
    if max_edges/min_edges > 100:
        max_edges = min_edges*100
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


def build_flow_graph(otherG, scaled_degree_in, scaled_degree_out):
    F = nx.DiGraph()

    for node in otherG.nodes():  # same nodes in pos/neg
        F.add_edge("source", ("out", node), capacity=scaled_degree_out.get(node, Fraction(0)).numerator)
        F.add_edge(("in", node), "sink", capacity=scaled_degree_in.get(node, Fraction(0)).numerator)

    for bait, prey in otherG.edges():
        F.add_edge(("out", bait), ("in", prey), capacity=1)

    return F


def write_degree_differance(filename, targetG, otherG, alpha_hat):
    in_degree, out_degree = get_degrees(targetG)
    in_degree_other, out_degree_other = get_degrees(otherG)

    with open(filename, "w") as w:
        w.write("node\ttarget_in_degree\tscaled_in_degree\ttarget_out_degree\tscaled_out_degree\testimated_increase\n")
        for node in targetG.nodes():
            w.write(f"{node}\t{in_degree.get(node, 0)}\t{in_degree_other.get(node, 0)}\t")
            w.write(f"{out_degree.get(node, 0)}\t{out_degree_other.get(node, 0)}\t{alpha_hat}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--positive_data", required=True, help="")
    parser.add_argument("--negative_data", required=True, help="")
    parser.add_argument("--max_flow_positive", required=True, help="Path to output csv file")
    parser.add_argument("--max_flow_negative", required=True, help="Path to output csv file")
    parser.add_argument("--min_max_flow", type=int, default=40, help="")
    parser.add_argument("--subset", type=str)

    args = parser.parse_args()
    positive_data = args.positive_data
    negative_data = args.negative_data
    max_flow_positive = args.max_flow_positive
    max_flow_negative = args.max_flow_negative
    min_max_flow = args.min_max_flow

    positive_bait_prey_df = pd.read_csv(positive_data, sep="\t")
    negative_bait_prey_df = pd.read_csv(negative_data, sep="\t")

    positive_bait_prey_df = positive_bait_prey_df[["gene_name_bait", "gene_name_prey"]]
    negative_bait_prey_df = negative_bait_prey_df[["gene_name_bait", "gene_name_prey"]]

    positive_bait_prey_df.columns = ["bait", "prey"]
    negative_bait_prey_df.columns = ["bait", "prey"]

    # baits = set(negative_bait_prey_df["bait"]) & set(positive_bait_prey_df["bait"])
    # all_prey = set(negative_bait_prey_df["prey"]) & set(positive_bait_prey_df["prey"])
    #
    # negative_bp_df = negative_bait_prey_df[
    #     negative_bait_prey_df["bait"].isin(baits) & negative_bait_prey_df["prey"].isin(all_prey)].copy()
    # positive_bp_df = positive_bait_prey_df[
    #     positive_bait_prey_df["bait"].isin(baits) & positive_bait_prey_df["prey"].isin(all_prey)].copy()

    positive_diG = nx.from_pandas_edgelist(
        positive_bait_prey_df, "bait", "prey", create_using=nx.DiGraph()
    )
    negative_diG = nx.from_pandas_edgelist(
        negative_bait_prey_df, "bait", "prey", create_using=nx.DiGraph()
    )
    success = False
    if args.subset == "test":
        pa = [Fraction(i, 1) for i in range(3,0,-1)]
    else:
        pa = [Fraction(i, 1) for i in range(10,0,-1)]
    for pai in pa:
        target_in, target_out, success = get_possible_scaling_factors(positive_diG, pai)
        if success:
            print(f"Trying a subset where {pai} : 1")
            testF = build_flow_graph(negative_diG, target_in, target_out)
            flow_value, flow_dict = nx.maximum_flow(testF, "source", "sink")
            percent_output = round(flow_value / sum(target_in.values()).numerator * 100)

            print(f"Flow value: {flow_value}, that being {percent_output} % of scaled degree")

            min_target_ppis = sum(target_in.values()) / pai
            min_ppi_target = min_target_ppis*.9 < flow_value < min_target_ppis*1.1

            if percent_output > min_max_flow or pai == 1 or min_ppi_target: # Fix this one later
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