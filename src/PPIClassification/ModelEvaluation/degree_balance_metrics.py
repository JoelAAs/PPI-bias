import argparse
import re
import pandas as pd
import numpy as np
import networkx as nx

def get_scaling_value(neg_file):
    with open(neg_file, "r") as f:
        first_line = f.readline()
    
    match = re.search(r"Scaled:\s*(\d+)", first_line)
    if match:
        return int(match.group(1))
    else:
        return None

def get_degree_df(degree_dict_pos, degree_dict_neg):
    df_d_pos = pd.DataFrame(degree_dict_pos, columns=["gene_name", "pos_degree"])
    df_d_neg = pd.DataFrame(degree_dict_neg, columns=["gene_name", "neg_degree"])
    
    degree_df=df_d_pos.merge(df_d_neg, on="gene_name", how="outer").fillna(0)

    return degree_df

def read_data_pair(pos_file, neg_file, network_type):
    if pos_file[-3:] == "csv":
        df_pos = pd.read_csv(pos_file, sep="\t", usecols=[0, 1], header=None)
        df_neg = pd.read_csv(neg_file, sep="\t", usecols=[0, 1], header=None, comment="#")        
        scale = get_scaling_value(neg_file)
        if not scale:
            scale=1
    else:
        df_pos = pd.read_parquet(pos_file).iloc[:,0:2]
        df_neg = pd.read_parquet(neg_file).iloc[:,0:2]
        scale=1

    graph_type = (nx.DiGraph if network_type == "directional" else nx.Graph)
    bp_cols = df_pos.columns
    G_pos = nx.from_pandas_edgelist(
        df_pos, bp_cols[0], bp_cols[1], create_using=graph_type
    )
    G_neg = nx.from_pandas_edgelist(
        df_neg, bp_cols[0], bp_cols[1], create_using=graph_type
    )

    return G_pos, G_neg, scale


def get_correlation(degree_df):
    pear_corr = degree_df["pos_degree"].corr(degree_df["neg_degree"], method="pearson")
    spear_corr = degree_df["pos_degree"].corr(degree_df["neg_degree"], method="spearman")

    return pear_corr, spear_corr

def get_wasserstein(degree_df, aimed_scaling):
    pos_scaled = (degree_df["pos_degree"] * aimed_scaling).astype(int)
    neg_degree = degree_df["neg_degree"].astype(int)

    max_degrees = int(max(pos_scaled.max(), neg_degree.max()))

    pos_hist = np.zeros(max_degrees)
    neg_hist = np.zeros(max_degrees)

    pos_counts = pos_scaled.value_counts()
    neg_counts = neg_degree.value_counts()

    pos_hist[pos_counts.index - 1] = pos_counts.values
    neg_hist[neg_counts.index - 1] = neg_counts.values

    pos_hist = pos_hist / pos_hist.sum()
    neg_hist = neg_hist / neg_hist.sum()
    cdf_pos = np.cumsum(pos_hist)
    cdf_neg = np.cumsum(neg_hist)

    wasserstein_distance = np.sum(np.abs(cdf_pos - cdf_neg))

    return wasserstein_distance

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get train test metrics")
    parser.add_argument("--pos_train", type=str, required=True)
    parser.add_argument("--neg_train", type=str, required=True)
    parser.add_argument("--pos_val", type=str, required=True)
    parser.add_argument("--neg_val", type=str, required=True)
    parser.add_argument("--pos_test", type=str, required=True)
    parser.add_argument("--neg_test", type=str, required=True)
    parser.add_argument("--output_file", type=str, required=True)
    parser.add_argument("--network_type", type=str, required=True)

    args = parser.parse_args()

    network_type=args.network_type
    pos_neg_pairs = zip ([
        args.pos_train,
        args.pos_val,
        args.pos_test
    ],[
        args.neg_train,
        args.neg_val,
        args.neg_test,
    ])
    with open(args.output_file, "w") as w:
        if network_type == "directional":
            w.write("dataset\tp_bait\tsp_bait\tws_bait\tp_prey\tsp_prey\tws_prey\taimed_scale\n")
            for p_f, n_f in pos_neg_pairs:
                G_pos, G_neg, aimed_scale = read_data_pair(p_f, n_f, network_type)
                
                filename = n_f.split("/")[-1]
                dataset = re.sub(r"_neg..*$", "", filename)
                
                degree_pos_bait = G_pos.out_degree()
                degree_neg_bait = G_neg.out_degree()
                degree_pos_prey = G_pos.in_degree()
                degree_neg_prey = G_neg.in_degree()
                    
                degree_bait_df = get_degree_df(degree_pos_bait, degree_neg_bait)
                degree_prey_df = get_degree_df(degree_pos_prey, degree_neg_prey)
                
                p_bait, sp_bait =  get_correlation(degree_bait_df)
                ws_bait = get_wasserstein(degree_bait_df, aimed_scale)
                p_prey, sp_prey =  get_correlation(degree_prey_df)
                ws_prey = get_wasserstein(degree_prey_df, aimed_scale)
                
                w.write(f"{dataset}\t{p_bait}\t{sp_bait}\t{ws_bait}\t{p_prey}\t{sp_prey}\t{ws_prey}\t{aimed_scale}\n")
        else:
            w.write("dataset\tp_undir\tsp_undir\tws_undir\taimed_scale\n")
            for p_f, n_f in pos_neg_pairs:
                G_pos, G_neg, aimed_scale = read_data_pair(p_f, n_f, network_type)
                filename = n_f.split("/")[-1]
                dataset = re.sub(r"_neg..*$", "", filename)
                
                degree_pos_undir = G_pos.degree()
                degree_neg_undir = G_neg.degree()
                
                degree_undir_df = get_degree_df(degree_pos_undir, degree_neg_undir)
                
                p_undir, sp_undir =  get_correlation(degree_undir_df)
                ws_undir = get_wasserstein(degree_undir_df, aimed_scale)
                
                w.write(f"{dataset}\t{p_undir}\t{sp_undir}\t{ws_undir}\t{aimed_scale}\n")
                
