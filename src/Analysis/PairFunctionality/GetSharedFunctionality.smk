import datetime

import numpy as np
import pandas as pd
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
    ji_index = np.zeros(3, float)
    intersects = np.zeros(3, int)
    for i, current_suffix in enumerate(term_suffix):
        intersect = int.bit_count(
            gene_idx_dicts[gene1][f"index_{current_suffix}"] & gene_idx_dicts[gene2][f"index_{current_suffix}"])
        union = int.bit_count(
            gene_idx_dicts[gene1][f"index_{current_suffix}"] | gene_idx_dicts[gene2][f"index_{current_suffix}"])
        ji_index[i] = intersect / union if union != 0 else np.nan
        intersects[i] = intersect

    return [*ji_index, *intersects]


rule get_jaccard_go_bait_prey:
    input:
        flat_df="work_folder/inferred_search_space/analysis/bias_reduced_ppis/p_estimated_protein_pairs.csv",
        abundance_df="work_folder/analysis/abundance_aware/POD_abundance.csv"
    output:
        flat_jaccard = "work_folder/analysis/GO/flat_jaccard.csv",
        abundance_jaccard = "work_folder/analysis/GO/abundance_jaccard.csv"

    run:
        go_terms = ["bp", "cc", "mf"]
        go_cols = [f"ji_{go}" for go in go_terms] + [f"intersect_{go}" for go in go_terms]

        flat_df = pd.read_csv(input.flat_df,sep="\t")
        all_flat_genes = set(flat_df["gene_name_bait"].unique()) | set(flat_df["gene_name_prey"].unique())
        gene_idx_dicts_flat = get_go_idx_df(all_flat_genes)
        flat_df[go_cols] = flat_df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard(x.iloc[0],x.iloc[1],gene_idx_dicts_flat),
            axis = 1,
            result_type='expand'
        )
        flat_df.to_csv(output.abundance_jaccard, sep="\t", index=False)
        del flat_df

        abundance_df = pd.read_csv(input.abundance_df,sep="\t")
        all_abundance_genes = list(
            set(abundance_df["gene_name_bait"].unique()) | set(abundance_df["gene_name_prey"].unique()))
        gene_idx_dicts_abundance = get_go_idx_df(all_abundance_genes)

        abundance_df[go_cols] = abundance_df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard(x.iloc[0],x.iloc[1],gene_idx_dicts_abundance),
            axis = 1,
            result_type='expand'
        )
        abundance_df.to_csv(output.abundance_jaccard, sep="\t", index=False)