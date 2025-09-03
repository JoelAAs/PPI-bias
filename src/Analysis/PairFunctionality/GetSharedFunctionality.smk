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



def get_sliding_avg_enrichment(df, value_column, greater=True, min_samples=50):
    bins = df[value_column].unique()
    bins.sort()
    values = df[value_column].values
    idx_val = values.argsort()
    if greater:
        idx_val = idx_val[::-1]
        bins = bins[::-1]
    values = values[idx_val]

    cc_ji = df["ji_bp"].to_numpy()[idx_val]
    bp_ji = df["ji_cc"].to_numpy()[idx_val]
    mf_ji = df["ji_mf"].to_numpy()[idx_val]
    cc_intersect = df["intersect_bp"].to_numpy()[idx_val]
    bp_intersect = df["intersect_cc"].to_numpy()[idx_val]
    mf_intersect = df["intersect_mf"].to_numpy()[idx_val]

    if not greater:
        bins = -bins
        values = -values

    cc_ji_sum = 0
    bp_ji_sum = 0
    mf_ji_sum = 0

    cc_ji_na = 0
    bp_ji_na = 0
    mf_ji_na = 0

    cc_intersect_sum = 0
    bp_intersect_sum = 0
    mf_intersect_sum = 0

    previous = 0
    i = 0
    j = 0
    rows = [{}] * len(bins)
    for threshold in bins:
        while (i < len(values) and threshold <= values[i]) or i < min_samples:
            i += 1

        cc_na_idx = np.isnan(cc_ji[previous:i])
        bp_na_idx = np.isnan(bp_ji[previous:i])
        mf_na_idx = np.isnan(mf_ji[previous:i])

        cc_ji_na += sum(cc_na_idx)
        bp_ji_na += sum(bp_na_idx)
        mf_ji_na += sum(mf_na_idx)

        cc_ji_sum += cc_ji[previous:i][~cc_na_idx].sum()
        bp_ji_sum += bp_ji[previous:i][~bp_na_idx].sum()
        mf_ji_sum += mf_ji[previous:i][~mf_na_idx].sum()
        cc_intersect_sum += cc_intersect[previous:i].sum()
        bp_intersect_sum += bp_intersect[previous:i].sum()
        mf_intersect_sum += mf_intersect[previous:i].sum()

        rows[j] = {
            "limit": value_column,
            "value": (threshold if greater else -threshold),
            "cc_ji_avg": cc_ji_sum / (i - cc_ji_na),
            "bp_ji_avg": bp_ji_sum / (i - bp_ji_na),
            "mf_ji_avg": mf_ji_sum / (i - mf_ji_na),
            "cc_intersect_avg": cc_intersect_sum / i,
            "bp_intersect_avg": bp_intersect_sum / i,
            "mf_intersect_avg": mf_intersect_sum / i,
            "number_of_pairs": i
        }
        if previous != i:
            j += 1
        previous = i

    return pd.DataFrame([r for r in rows if r])

def get_pair_jaccard(gene1, gene2, gene_idx_dicts):
    term_suffix = ["bp", "cc", "mf"]
    ji_index = np.zeros(3,float)
    intersects = np.zeros(3,int)
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
        flat_jaccard="work_folder/analysis/GO/flat_jaccard.csv",
        abundance_jaccard="work_folder/analysis/GO/abundance_jaccard.csv"

    run:
        go_terms = ["bp", "cc", "mf"]
        go_cols = [f"ji_{go}" for go in go_terms] + [f"intersect_{go}" for go in go_terms]

        flat_df = pd.read_csv(input.flat_df,sep="\t")
        all_flat_genes = set(flat_df["gene_name_bait"].unique()) | set(flat_df["gene_name_prey"].unique())
        gene_idx_dicts_flat = get_go_idx_df(all_flat_genes)
        flat_df[go_cols] = flat_df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard(x.iloc[0],x.iloc[1],gene_idx_dicts_flat),
            axis=1,
            result_type='expand'
        )
        flat_df.to_csv(output.flat_jaccard,sep="\t",index=False)
        del flat_df

        abundance_df = pd.read_csv(input.abundance_df,sep="\t")
        all_abundance_genes = list(
            set(abundance_df["gene_name_bait"].unique()) | set(abundance_df["gene_name_prey"].unique()))
        gene_idx_dicts_abundance = get_go_idx_df(all_abundance_genes)

        abundance_df[go_cols] = abundance_df[["gene_name_bait", "gene_name_prey"]].apply(
            lambda x: get_pair_jaccard(x.iloc[0],x.iloc[1],gene_idx_dicts_abundance),
            axis=1,
            result_type='expand'
        )
        abundance_df.to_csv(output.abundance_jaccard,sep="\t",index=False)


rule get_go_accumulation:
    input:
        flat_jaccard="work_folder/analysis/GO/flat_jaccard.csv",
        abundance_jaccard="work_folder/analysis/GO/abundance_jaccard.csv"
    output:
        flat_jaccard_greater="work_folder/analysis/GO/flat_jaccard_greater.csv",
        abundance_jaccard_greater="work_folder/analysis/GO/abundance_jaccard_greater.csv",
        flat_jaccard_lesser="work_folder/analysis/GO/flat_jaccard_lesser.csv",
        abundance_jaccard_lesser="work_folder/analysis/GO/abundance_jaccard_lesser.csv"
    run:
        abundance_df = pd.read_csv(
            input.abundance_jaccard, sep="\t"
        )
        get_sliding_avg_enrichment(
            abundance_df,
            "lower_bound_pod").to_csv(
            output.abundance_jaccard_greater,
            sep="\t",index=False
        )
        get_sliding_avg_enrichment(
            abundance_df,
            "upper_bound_pod",
            greater=False).to_csv(
            output.abundance_jaccard_lesser,
            sep="\t", index=False
        )

        flat_df = pd.read_csv(
            input.flat_jaccard, sep="\t"
        )
        get_sliding_avg_enrichment(
            flat_df,
            "p_lower_ci").to_csv(
            output.flat_jaccard_greater,
            sep="\t",index=False
        )
        get_sliding_avg_enrichment(
            flat_df,
            "p_upper_ci",
            greater=False).to_csv(
            output.flat_jaccard_lesser,
            sep="\t", index=False
        )


rule plot_go_accumulation:
    input:
        flat_jaccard_greater = "work_folder/analysis/GO/flat_jaccard_greater.csv",
        abundance_jaccard_greater = "work_folder/analysis/GO/abundance_jaccard_greater.csv",
        flat_jaccard_lesser = "work_folder/analysis/GO/flat_jaccard_lesser.csv",
        abundance_jaccard_lesser = "work_folder/analysis/GO/abundance_jaccard_lesser.csv"
    output:
        go_jaccard = "work_folder/plots/GO/intersect_GO_vs_POD.png",
        go_intersect = "work_folder/plots/GO/jaccard_GO_vs_POD.png"
    shell:
        """
        Rscript src/Plotting/plot_go_accumulation.R
        """
