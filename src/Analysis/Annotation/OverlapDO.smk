import numpy as np
import pandas as pd
from mean_distance_support import get_cumulative_sum
from mygene import MyGeneInfo


def get_mask(obs_do, all_dos):
    mask = np.zeros(len(all_dos),dtype=bool)
    i = 0
    while i < len(all_dos) and len(obs_do) > 0:
        mask[i] = all_dos[i] in obs_do
        if mask[i]:
            obs_do.remove(all_dos[i])
        i += 1
    return mask


def get_base10_idx(mask):
    mask = mask.astype(int).astype(str)
    res = int("".join(mask),2)
    return res


def get_do_idx_df(do_file):
    do_df = pd.read_csv(do_file,sep="\t")
    do_per_gene_df = do_df[["gene_name", "doid"]].groupby("gene_name",as_index=False)["doid"].unique()
    all_dos = do_df["doid"].dropna().unique()
    do_per_gene_dict = {gene: doid[doid == doid].tolist() for i, (gene, doid) in do_per_gene_df.iterrows()}

    gene_idx_dicts = dict()
    for gene in do_per_gene_dict.keys():
        do_idx = get_base10_idx(get_mask(do_per_gene_dict[gene].copy(),all_dos))
        gene_idx_dicts[gene] = do_idx
    return gene_idx_dicts


def get_pair_jaccard_do(gene1, gene2, gene_idx_dicts):
    intersect = int.bit_count(
        gene_idx_dicts[gene1] & gene_idx_dicts[gene2])
    union = int.bit_count(
        gene_idx_dicts[gene1] | gene_idx_dicts[gene2])
    ji_index = intersect / union if union != 0 else np.nan
    n_do_prey = int.bit_count(gene_idx_dicts[gene1])
    n_do_bait = int.bit_count(gene_idx_dicts[gene2])

    return [ji_index, intersect, n_do_prey, n_do_bait]


rule get_gene_do_terms:
    params:
        script="src/Analysis/Annotation/get_DO.R"
    input:
        pod_df=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        do_terms=f"work_folder{pn}/analysis/DO/gene_do_{{data}}.txt"
    conda: "do_enrichment"
    shell:
        """
        Rscript {params.script} {input.pod_df} {output.do_terms}
        """

rule get_jaccard_do_bait_prey:
    input:
        do_gene_df=f"work_folder{pn}/analysis/DO/gene_do_{{data}}.txt",
        pod_df=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
    output:
        do_jaccard=f"work_folder{pn}/analysis/DO/pairs_{{data}}_jaccard.csv"
    run:
        do_cols = [
            "ji_do",
            "intersect_do",
            "n_do_bait",
            "n_do_prey"
        ]

        df = pd.read_csv(input.pod_df,sep="\t")
        all_genes = set(df["gene_name_bait"].unique()) | set(df["gene_name_prey"].unique())
        gene_idx_dicts_flat = get_do_idx_df(input.do_gene_df)
        df[do_cols] = df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard_do(x.iloc[0],x.iloc[1],gene_idx_dicts_flat),
            axis=1,
            result_type='expand'
        )
        df[["pair_id"] + do_cols].to_csv(output.do_jaccard,sep="\t",index=False)


rule get_do_accumulation:
    input:
        pod_df=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
        do_jaccard=f"work_folder{pn}/analysis/DO/pairs_{{data}}_jaccard.csv"
    output:
        jaccard_greater=f"work_folder{pn}/analysis/DO/cumulative/POD_{{data}}_jaccard_greater.csv",
        jaccard_lesser=f"work_folder{pn}/analysis/DO/cumulative/POD_{{data}}_jaccard_lesser.csv"
    run:
        do_data = pd.read_csv(
            input.do_jaccard,sep="\t"
        )
        pod_df = pd.read_csv(
            input.pod_df,
            sep="\t"
        )
        do_df = pod_df.merge(do_data, on="pair_id")

        do_cols = [
            "ji_do",
            "intersect_do",
            "n_do_bait",
            "n_do_prey"
        ]

        get_cumulative_sum(
            do_df,
            value_column="lower_bound_pod",
            cumulative_columns=do_cols
        ).to_csv(
            output.jaccard_greater,
            sep="\t",index=False
        )
        get_cumulative_sum(
            do_df,
            value_column="upper_bound_pod",
            cumulative_columns=do_cols,
            greater=False
        ).to_csv(
            output.jaccard_lesser,
            sep="\t",index=False
        )
