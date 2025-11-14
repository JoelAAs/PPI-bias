import numpy as np
import pandas as pd
from mean_distance_support import get_cumulative_sum
from mygene import MyGeneInfo


def get_go_genes(genes):
    mg = MyGeneInfo()
    result = mg.querymany(
        genes,
        scopes="symbol",
        fields="go",
        species="human",
        returnall=True
    )

    gene_go_dict = dict()
    for go_q in result["out"]:
        gene_gos = dict()
        gene = go_q["query"]
        go_terms = go_q.get('go',{})
        for i, term in enumerate(["BP", "MF", "CC"]):
            go_match = go_terms.get(term,{})
            if isinstance(go_match,list):
                gos = {go["id"] for go in go_match}
            elif "id" in go_match:
                gos = {go_match["id"], }
            else:
                gos = set()
            gene_gos[term] = gos

        gene_go_dict[gene] = gene_gos

    for missing_gene in result["missing"]:
        gene_gos = {
            "BP": set(),
            "CC": set(),
            "MF": set()
        }
        gene_go_dict[missing_gene] = gene_gos
    return gene_go_dict


def get_mask(obs_go, all_gos):
    mask = np.zeros(len(all_gos),dtype=bool)
    i = 0
    while i < len(all_gos) and obs_go:
        mask[i] = all_gos[i] in obs_go
        if mask[i]:
            obs_go.remove(all_gos[i])
        i += 1
    return mask


def get_base10_idx(mask):
    mask = mask.astype(int).astype(str)
    res = int("".join(mask),2)
    return res


def get_go_idx_df(genes):
    go_per_gene = get_go_genes(genes)
    all_bp = sorted(list(set.union(*[go_per_gene[gene]["BP"] for gene in go_per_gene.keys()])))
    all_cc = sorted(list(set.union(*[go_per_gene[gene]["CC"] for gene in go_per_gene.keys()])))
    all_mf = sorted(list(set.union(*[go_per_gene[gene]["MF"] for gene in go_per_gene.keys()])))

    gene_idx_dicts = dict()
    for gene in go_per_gene.keys():
        bp_idx = get_base10_idx(get_mask(go_per_gene[gene]["BP"].copy(),all_bp))
        cc_idx = get_base10_idx(get_mask(go_per_gene[gene]["CC"].copy(),all_cc))
        mf_idx = get_base10_idx(get_mask(go_per_gene[gene]["MF"].copy(),all_mf))

        gene_idx_dicts[gene] = {
            "index_bp": bp_idx,
            "index_cc": cc_idx,
            "index_mf": mf_idx
        }
    return gene_idx_dicts


def get_pair_jaccard(gene1, gene2, gene_idx_dicts):
    term_suffix = ["bp", "cc", "mf"]
    ji_index = np.zeros(3,float)
    intersects = np.zeros(3,int)
    n_go_prey = np.zeros(3,float)
    n_go_bait = np.zeros(3,float)

    for i, current_suffix in enumerate(term_suffix):
        intersect = int.bit_count(
            gene_idx_dicts[gene1][f"index_{current_suffix}"] & gene_idx_dicts[gene2][f"index_{current_suffix}"])
        union = int.bit_count(
            gene_idx_dicts[gene1][f"index_{current_suffix}"] | gene_idx_dicts[gene2][f"index_{current_suffix}"])
        ji_index[i] = intersect / union if union != 0 else np.nan
        intersects[i] = intersect
        n_go_prey[i] = int.bit_count(gene_idx_dicts[gene1][f"index_{current_suffix}"])
        n_go_bait[i] = int.bit_count(gene_idx_dicts[gene2][f"index_{current_suffix}"])

    return [*ji_index, *intersects, *n_go_prey, *n_go_bait]


rule get_jaccard_go_bait_prey:
    input:
        pod_df=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        go_jaccard=f"work_folder{pn}/analysis/GO/pairs_{{data}}_jaccard.csv"
    run:
        go_terms = ["bp", "cc", "mf"]
        go_cols = (
                [f"ji_{go}" for go in go_terms] +
                [f"intersect_{go}" for go in go_terms] +
                [f"n_go_bait_{go}" for go in go_terms] +
                [f"n_go_prey_{go}" for go in go_terms]
        )

        df = pd.read_csv(input.pod_df,sep="\t")
        all_genes = set(df["gene_name_bait"].unique()) | set(df["gene_name_prey"].unique())
        gene_idx_dicts_flat = get_go_idx_df(all_genes)
        df[go_cols] = df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard(x.iloc[0],x.iloc[1],gene_idx_dicts_flat),
            axis=1,
            result_type='expand'
        )
        df[["pair_id"] + go_cols].to_csv(output.go_jaccard,sep="\t",index=False)


rule get_go_accumulation:
    input:
        go_jaccard=f"work_folder{pn}/analysis/GO/pairs_{{data}}_jaccard.csv",
        pod_df=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        jaccard_greater=f"work_folder{pn}/analysis/GO/cumulative/POD_{{data}}_jaccard_greater.csv",
        jaccard_lesser=f"work_folder{pn}/analysis/GO/cumulative/POD_{{data}}_jaccard_lesser.csv"
    run:
        go_data = pd.read_csv(
            input.go_jaccard,sep="\t"
        )

        pod_df = pd.read_csv(
            input.pod_df,
            sep="\t"
        )
        go_df = pod_df.merge(go_data,on="pair_id")

        go_terms = ["bp", "cc", "mf"]
        go_cols = (
                [f"ji_{go}" for go in go_terms] +
                [f"intersect_{go}" for go in go_terms] +
                [f"n_go_bait_{go}" for go in go_terms] +
                [f"n_go_prey_{go}" for go in go_terms]
        )


        get_cumulative_sum(
            go_df,
            value_column="lower_bound_pod",
            cumulative_columns=go_cols
        ).to_csv(
            output.jaccard_greater,
            sep="\t",index=False
        )
        get_cumulative_sum(
            go_df,
            value_column="upper_bound_pod",
            cumulative_columns=go_cols,
            greater=False
        ).to_csv(
            output.jaccard_lesser,
            sep="\t",index=False
        )


rule abundance_go_plot:
    input:
        norm_log="data/normalised_log_ra.csv"
    output:
        goterms_abundance=f"work_folder{pn}/analysis/GO/ra_pod_vs_go_terms.csv"
    run:
        abundance_df = pd.read_csv(input.norm_log,sep="\t")
        del abundance_df["samples"]
        del abundance_df["cell_line"]

        abundance_mean = abundance_df.mean(axis=0,skipna=True).reset_index()
        abundance_mean.columns = ["gene_name", "relative_abundance"]
        abundance_mean["pod"] = (~abundance_df.isna()).mean(axis=0).values
        go_dict = get_go_genes(abundance_mean["gene_name"])


        def get_n_go(gene, go_dict):
            return [
                len(go_dict[gene]["BP"]),
                len(go_dict[gene]["MF"]),
                len(go_dict[gene]["CC"])
            ]


        abundance_mean[[
            "n_bp",
            "n_mf",
            "n_cc"
        ]] = abundance_mean.apply(
            lambda x: get_n_go(x["gene_name"],go_dict),axis=1,result_type='expand')

        abundance_mean.to_csv(output.goterms_abundance,sep="\t",index=False)


rule bait_usage:
    input:
        bait_prey=f"work_folder{pn}/formated/bait_prey_publications.csv"
    output:
        goterms_studies=f"work_folder{pn}/analysis/GO/n_studies_go_terms.csv"
    run:
        df = pd.read_csv(input.bait_prey,sep="\t")
        df_bait = df[
            ~df[["pubmed_id", "detection_method", "gene_name_bait"]].duplicated()
        ].groupby("gene_name_bait",as_index=False).size()
        df_bait.columns = ["gene_name", "n_studies"]
        go_dict = get_go_genes(df_bait["gene_name"])


        def get_n_go(gene, go_dict):
            return [
                len(go_dict[gene]["BP"]),
                len(go_dict[gene]["MF"]),
                len(go_dict[gene]["CC"])
            ]


        df_bait[[
            "n_bp",
            "n_mf",
            "n_cc"
        ]] = df_bait.apply(
            lambda x: get_n_go(x["gene_name"],go_dict),axis=1,result_type='expand')
        df_bait.to_csv(output.goterms_studies,sep="\t",index=False)


