import pandas as pd
import numpy as np


def get_nrow(file):
    with open(file, "r") as f:
        nlines = sum(1 for _ in f)
    return nlines


def pos_neg_degree_cor(posdf, negdf, degree_name):
    pos_degree = posdf.groupby(degree_name, as_index=False).size()
    neg_degree = negdf.groupby(degree_name, as_index=False).size()
    deg_df = pd.merge(pos_degree, neg_degree, on=degree_name, how="outer", suffixes=("_p", "_n")).fillna(0)
    deg_df[["size_p", "size_n"]].corr(method="pearson")
    pearson = deg_df[["size_p", "size_n"]].corr(method="pearson").iloc[0, 1]
    spearman = deg_df[["size_p", "size_n"]].corr(method="spearman").iloc[0, 1]
    return spearman, pearson


def get_degree_difference(file_pattern):
    if file_pattern[-3] == "csv":
        read_file = lambda x: pd.read_csv(x, sep="\t")
    else:
        read_file = lambda x: pd.read_csv(x, sep="\t", comment="#", header=None).rename({0: "bait", 1: "prey"}, axis=1)

    pos_df = read_file(file_pattern.format(datatype="pos"))
    neg_df = read_file(file_pattern.format(datatype="neg"))
    n_bait_pos = pos_df.bait.nunique()
    n_prey_pos = pos_df.prey.nunique()
    n_bait_neg = neg_df.bait.nunique()
    n_prey_neg = neg_df.prey.nunique()

    sb, pb = pos_neg_degree_cor(pos_df, neg_df, "bait")
    sp, pb = pos_neg_degree_cor(pos_df, neg_df, "prey")

    return [
        ["bait", sb, pb, n_bait_pos, n_bait_neg],
        ["prey", sp, sp, n_prey_pos, n_prey_neg]
    ]


def get_order(dataset, typeset, datatype, pn="/per_gene"):
    order = (
        f"work_folder{pn}/subsets/{dataset}_full_{datatype}.csv",
        f"work_folder{pn}/subsets/{typeset}/{dataset}_{datatype}.csv",
        f"work_folder{pn}/subsets/{typeset}/maxflow/{dataset}_{datatype}.edgelist",
        f"work_folder{pn}/subsets/{typeset}/balanced/{dataset}_{datatype}.csv"
    )
    n_rows = [get_nrow(o) for o in order]
    return n_rows


def write_number_of_edges(subsetting_file, dataset):
    with open(subsetting_file, "w") as w:
        for typeset in ["train", "test", "validation"]:
            for datatype in ["pos", "neg"]:
                line = get_order(dataset, typeset, datatype, pn="/per_gene")
                w.write("\t".join([dataset, typeset, datatype]) + "\t" +
                        "\t".join(map(str, line)) + "\n")


def get_cor_files(dataset, cor_file, pn="/per_gene"):
    datatype = "{datatype}"
    with open(cor_file, "w") as w:
        # full
        lines = get_degree_difference(
            f"work_folder{pn}/subsets/{dataset}_full_{datatype}.csv",
        )
        for line in lines:
            w.write(f"Full\tAll\t" + "\t".join(map(str, line)) + "\n")

        for typeset in ["train", "test", "validation"]:
            # kaffpa
            lines = get_degree_difference(
                f"work_folder{pn}/subsets/{typeset}/{dataset}_{datatype}.csv"
            )
            for line in lines:
                w.write(f"kaffpa\t{typeset}\t" + "\t".join(map(str, line)) + "\n")
            # maxflow
            lines = get_degree_difference(
                f"work_folder{pn}/subsets/{typeset}/maxflow/{dataset}_{datatype}.edgelist"
            )
            for line in lines:
                w.write(f"maxflow\t{typeset}\t" + "\t".join(map(str, line)) + "\n")
            # ilp
            lines = get_degree_difference(
                f"work_folder{pn}/subsets/{typeset}/balanced/{dataset}_{datatype}.csv"
            )
            for line in lines:
                w.write(f"ILP\t{typeset}\t" + "\t".join(map(str, line)) + "\n")


if __name__ == '__main__':
    for subset in ["ms", "y2h", "flat"]:
        get_cor_files(subset, f"work_folder/per_gene/subsets/summary_statistics/correlation_{subset}.csv", pn="/per_gene")
        write_number_of_edges(f"work_folder/per_gene/subsets/summary_statistics/n_edges_{subset}.csv", subset)