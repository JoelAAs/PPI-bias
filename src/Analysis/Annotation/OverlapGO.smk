import numpy as np
import pandas as pd
from mean_distance_support import get_cumulative_sum
from mygene import MyGeneInfo
from collections import Counter


def get_go_genes(genes):
    mg = MyGeneInfo()
    result = mg.querymany(
        genes,
        scopes="symbol,alias",
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
        pod_df="work_folder/analysis/POD/POD_{data}.csv"
    output:
        go_jaccard="work_folder/analysis/GO/pairs_{data}_jaccard.csv"
    log:
        "logs/analysis/GO/pairs_{data}_jaccard.log"
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
        go_jaccard="work_folder/analysis/GO/pairs_{data}_jaccard.csv",
        pod_df="work_folder/analysis/POD/POD_{data}.csv"
    output:
        jaccard_greater="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_greater.csv",
        jaccard_lesser="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_lesser.csv"
    log:
        "logs/analysis/GO/cumulative/POD_{data}_jaccard.log"
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
        goterms_abundance="work_folder/analysis/GO/ra_pod_vs_go_terms.csv"
    log:
        "logs/analysis/GO/ra_pod_vs_go_terms.log"
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
        bait_prey="work_folder/formated/bait_prey_publications.csv"
    output:
        goterms_studies="work_folder/analysis/GO/n_studies_go_terms.csv"
    log:
        "logs/analysis/GO/n_studies_go_terms.log"
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


rule get_huri_mf_counts:
    input:
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063.csv",
        huri="work_folder/data/huri/intact_huri.csv"
    output:
        compare_data="work_folder/analysis/GO/huri_vs_bioplex_annotation.csv",
        summary_data="work_folder/analysis/GO/huri_vs_bioplex_shared_ji.csv"
    log:
        "logs/analysis/GO/huri_vs_bioplex_annotation.log"
    run:
        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]

        huri_df = pd.read_csv(input.huri,sep="\t")[["gene_name_bait", "gene_name_prey"]]
        #shared_baits = set(bp_df["gene_name_bait"]) & set(huri_df["gene_name_bait"])

        #huri_df = huri_df[huri_df["gene_name_bait"].isin(shared_baits)]
        #bp_df = bp_df[bp_df["gene_name_bait"].isin(shared_baits)]

        all_genes = set(bp_df["gene_name_bait"]) | \
                    set(bp_df["gene_name_prey"]) | \
                    set(huri_df["gene_name_bait"]) | \
                    set(huri_df["gene_name_prey"])

        go_dict = get_go_genes(all_genes)


        def flat_and_count(nested_list, method, role):
            n_ppis = len(nested_list)
            flat_list = [i for l in nested_list for i in l]
            count = Counter(flat_list)
            mf_df = pd.DataFrame(count.items(),columns=["mf_terms", "count"])
            mf_df["mf_frequency"] = mf_df["count"] / n_ppis
            mf_df[["method", "role"]] = [method, role]

            return mf_df


        def get_overlap_frequency(df, method):
            df["shared_mf"] = df.apply(lambda x: x["mf_terms_bait"] & x["mf_terms_prey"],axis=1)
            df["union_mf"] = df.apply(lambda x: x["mf_terms_bait"] | x["mf_terms_prey"],axis=1)
            flat_list = [i for l in df["shared_mf"].tolist() for i in l]
            count = Counter(flat_list)
            mf_df = pd.DataFrame(count.items(),columns=["mf_terms", "count"])
            mf_df["mf_frequency"] = mf_df["count"] / df.shape[0]
            mf_df[["method", "role"]] = [method, "Shared"]
            return mf_df, df


        huri_df["mf_terms_bait"] = huri_df["gene_name_bait"].apply(lambda x: go_dict[x]["MF"])
        huri_df["mf_terms_prey"] = huri_df["gene_name_prey"].apply(lambda x: go_dict[x]["MF"])
        bp_df["mf_terms_bait"] = bp_df["gene_name_bait"].apply(lambda x: go_dict[x]["MF"])
        bp_df["mf_terms_prey"] = bp_df["gene_name_prey"].apply(lambda x: go_dict[x]["MF"])

        huri_shared_mf_count, huri_df = get_overlap_frequency(huri_df,"HuRi")
        bp_shared_mf_count, bp_df = get_overlap_frequency(bp_df,"Bioplex")

        huri_bait_mf_count = flat_and_count(
            huri_df["mf_terms_bait"],
            "HuRi",
            "Bait"
        )
        huri_prey_mf_count = flat_and_count(
            huri_df["mf_terms_prey"],
            "HuRi",
            "Prey"
        )

        bp_bait_mf_count = flat_and_count(
            bp_df["mf_terms_bait"],
            "Bioplex",
            "Bait"
        )
        bp_prey_mf_count = flat_and_count(
            bp_df["mf_terms_prey"],
            "Bioplex",
            "Prey"
        )

        mf_frequency_df = pd.concat(
            [
                huri_shared_mf_count,
                bp_shared_mf_count,
                huri_bait_mf_count,
                huri_prey_mf_count,
                bp_bait_mf_count,
                bp_prey_mf_count
            ]
        )
        mf_frequency_df.to_csv(output.compare_data,sep="\t",index=None)


        def get_ji(row):
            if len(row["union_mf"]) != 0:
                return len(row["shared_mf"]) / len(row["union_mf"])
            else:
                return None


        huri_df["jaccard"] = huri_df.apply(get_ji,axis=1)
        bp_df["jaccard"] = bp_df.apply(get_ji,axis=1)
        without_mf_huri = sum((huri_df["mf_terms_bait"].apply(len) == 0) | (huri_df["mf_terms_prey"].apply(len) == 0))
        without_mf_bp = sum((bp_df["mf_terms_bait"].apply(len) == 0) | (bp_df["mf_terms_prey"].apply(len) == 0))

        with open(output.summary_data,"w") as w:
            w.write(f"HuRi jaccard index: {huri_df["jaccard"].mean()}\n")
            w.write(f"Bioplex jaccard index: {bp_df["jaccard"].mean()}\n")
            w.write(f"HuRi mean union: {huri_df["union_mf"].apply(len).mean()}\n")
            w.write(f"Bioplex mean union: {bp_df["union_mf"].apply(len).mean()}\n")
            w.write(f"HuRi mean intersect: {huri_df["shared_mf"].apply(len).mean()}\n")
            w.write(f"Bioplex mean intersect: {bp_df["shared_mf"].apply(len).mean()}\n")
            w.write(f"HuRi no annotation: {without_mf_huri}/{huri_df.shape[0]}\n")
            w.write(f"Bioplex no annotations: {without_mf_bp}/{bp_df.shape[0]}\n")
