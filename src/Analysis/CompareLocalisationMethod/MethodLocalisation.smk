import mygene
import re
from UniProtMapper import ProtMapper
import pandas as pd
import numpy as np
from src.Analysis.Annotation.localisation_support import add_localisation

rule get_huri:
    """
    must be annotated on gene name as bioplex is reported on gene name 
    """
    input:
        intact=f"work_folder{pn}/formated/bait_prey_publications.csv"
    output:
        huri=f"work_folder{pn}/data/huri/intact_huri.csv"
    run:
        intact = pd.read_csv(input.intact,sep="\t")
        huri_df = intact[intact["pubmed_id"] == 32296183]
        huri_df = huri_df[huri_df["gene_name_prey"] != huri_df["gene_name_bait"]] # Not removed if we are interested in isofroms
        huri_df.to_csv(output.huri,sep="\t",index=False)

rule localisation_delta:
    # TODO: reformat
    input:
        cvcl_0063_bp=f"work_folder{pn}/data/bioplex/CVCL_0063.csv",
        huri=f"work_folder{pn}/data/huri/intact_huri.csv",
        localisation_csv=f"work_folder{pn}/analysis/localisation/gene_to_localisation.csv"
    output:
        localisation_method=f"work_folder{pn}/analysis/localisation/HuRI_vs_Bioplex.csv"
    run:
        n_permutations = 100000

        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        huri_df = pd.read_csv(input.huri,sep="\t")
        df_localisation = pd.read_csv(input.localisation_csv,sep="\t")

        bp_df = add_localisation(bp_df,df_localisation)
        huri_df = add_localisation(huri_df,df_localisation)
        bait_intersect = set(bp_df["gene_name_bait"]) & set(huri_df["gene_name_bait"])
        shared_baits_bp_df = bp_df[bp_df["gene_name_bait"].isin(bait_intersect)]
        shared_baits_huri_df = huri_df[huri_df["gene_name_bait"].isin(bait_intersect)]

        # Prey localisation PPI count:
        local_bait_huri = shared_baits_huri_df.groupby(
            ["gene_name_bait", "localisation_bait"],as_index=False).size().rename({"size": "huri_ppi"},axis=1)
        local_bait_bp = shared_baits_bp_df.groupby(
            ["gene_name_bait", "localisation_bait"],as_index=False).size().rename({"size": "bp_ppi"},axis=1)

        local_bait_ppis = local_bait_bp.merge(
            local_bait_huri,on=["gene_name_bait", "localisation_bait"],how="outer").fillna(0)

        keep_local = local_bait_ppis.groupby("localisation_bait").size()
        keep_local = keep_local[keep_local > 10].index.values
        local_bait_ppis = local_bait_ppis[local_bait_ppis["localisation_bait"].isin(keep_local)]

        per_mat = np.zeros(shape=(n_permutations, len(keep_local)))
        sample_size = round((1 - 0.10) * local_bait_ppis.shape[0])
        for i in range(n_permutations):  # SD on baits not individual ppis
            ss_bait_ppis = local_bait_ppis.sample(sample_size)
            ss_bait_ppis["huri_ppi"] = ss_bait_ppis["huri_ppi"]/ss_bait_ppis["huri_ppi"].sum()
            ss_bait_ppis["bp_ppi"] = ss_bait_ppis["bp_ppi"]/ss_bait_ppis["bp_ppi"].sum()
            ss_bait_ppis["delta"] = ss_bait_ppis["huri_ppi"] - ss_bait_ppis["bp_ppi"]
            per_mat[i, :] = ss_bait_ppis.groupby("localisation_bait")["delta"].sum().values

        local_bait_ppis["huri_ppi"] = local_bait_ppis["huri_ppi"]/local_bait_ppis["huri_ppi"].sum()
        local_bait_ppis["bp_ppi"] = local_bait_ppis["bp_ppi"]/local_bait_ppis["bp_ppi"].sum()
        local_bait_ppis["delta"] = local_bait_ppis["huri_ppi"] - local_bait_ppis["bp_ppi"]
        bait_summed_delta = local_bait_ppis.groupby("localisation_bait",as_index=False)["delta"].sum()
        bait_summed_delta["localisation"] = bait_summed_delta["localisation_bait"]
        bait_summed_delta[["ci_2.5", "ci_97.5"]] = np.percentile(per_mat,[2.5, 97.5],axis=0).transpose()
        bait_summed_delta["role"] = "Bait"

        # PPI count non overlapping
        localisation_bp = shared_baits_bp_df.groupby(
            ["gene_name_prey", "localisation_prey"],as_index=False).size().rename({
            "size": "bp_ppi"},axis=1)
        localisation_huri = shared_baits_huri_df.groupby(
            ["gene_name_prey", "localisation_prey"],as_index=False).size().rename({
            "size": "huri_ppi"},axis=1)
        localisation_prey = localisation_bp.merge(localisation_huri,
            on=["gene_name_prey", "localisation_prey"],how="outer").fillna(0)

        keep_local = localisation_prey.groupby("localisation_prey").size()
        keep_local = keep_local[keep_local > 10].index.values
        localisation_prey = localisation_prey[localisation_prey["localisation_prey"].isin(keep_local)]

        per_mat = np.zeros(shape=(n_permutations, len(keep_local)))
        sample_size = round((1 - 0.10) * localisation_prey.shape[0])
        for i in range(n_permutations):  # SD on baits not individual ppis
            ss_prey_ppis = localisation_prey.sample(sample_size)
            ss_prey_ppis["huri_ppi"] = ss_prey_ppis["huri_ppi"]/ss_prey_ppis["huri_ppi"].sum()
            ss_prey_ppis["bp_ppi"] = ss_prey_ppis["bp_ppi"]/ss_prey_ppis["bp_ppi"].sum()
            ss_prey_ppis["delta"] = ss_prey_ppis["huri_ppi"] - ss_prey_ppis["bp_ppi"]
            prey_delta = ss_prey_ppis.groupby("localisation_prey")["delta"].sum()
            ordered = prey_delta.reindex(keep_local)
            per_mat[i, :] = ordered

        localisation_prey["huri_ppi"] = localisation_prey["huri_ppi"]/localisation_prey["huri_ppi"].sum()
        localisation_prey["bp_ppi"] = localisation_prey["bp_ppi"]/localisation_prey["bp_ppi"].sum()
        localisation_prey["delta"] = localisation_prey["huri_ppi"] - localisation_prey["bp_ppi"]
        prey_summed_delta = localisation_prey.groupby("localisation_prey",as_index=False)["delta"].sum()
        prey_summed_delta["localisation"] = prey_summed_delta["localisation_prey"]
        prey_summed_delta[["ci_2.5", "ci_97.5"]] = np.nanpercentile(per_mat,[2.5, 97.5],axis=0).transpose()
        prey_summed_delta["role"] = "Prey"
        columns = ['delta', 'ci_2.5', 'ci_97.5', 'role', 'localisation']
        pd.concat(
            [prey_summed_delta[columns], bait_summed_delta[columns]]
        ).to_csv(output.localisation_method,sep="\t",index=None)


def rec_nested_get(keys, dict):
    if len(keys) == 1:
        return dict.get(keys[0],[])
    elif keys[0] in dict:
        res = rec_nested_get(keys[1:],dict[keys[0]])
        if type(res) is not list:
            return [res, ]
        else:
            return res
    return []


def map_membrane_structure(gene_names):
    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        gene_names,
        scopes='symbol, alias',
        fields='uniprot.Swiss-Prot',
        species='human')

    rows = [
        (result["query"], prot_id)
        for result in results for prot_id in rec_nested_get(["uniprot", "Swiss-Prot"],result)
    ]
    map_df = pd.DataFrame(rows,columns=["gene_name", "uniprot_id"])
    map_df = map_df.dropna()
    return map_df


def get_membrane_annotation(map_df):
    mapper = ProtMapper()
    feature_columns = ["ft_region", "ft_transmem", "organism_name", "annotation_score"]
    result_df, failed = mapper.get(
        ids=map_df["uniprot_id"].unique(),
        from_db='UniProtKB_AC-ID',
        fields=feature_columns)

    contains_membrane = lambda x: bool(re.search(r"[Mm]embrane",x))
    result_df["membrane_association"] = result_df["Region"].apply(contains_membrane) | result_df[
        "Transmembrane"].apply(bool)
    trans_df = map_df.merge(result_df,left_on="uniprot_id",right_on="From")
    trans_df = trans_df.sort_values("Annotation",ascending=False)
    trans_df = trans_df.drop_duplicates(subset="gene_name",keep="first")  # keep the first with the highest annotation

    return trans_df[["gene_name", "membrane_association"]]


def join_membrane_annotation(ppi_df, membrane_df):
    ppi_df = ppi_df.merge(
        membrane_df,
        left_on="gene_name_bait",
        right_on="gene_name",
        how="left")
    del ppi_df["gene_name"]

    ppi_df = ppi_df.merge(
        membrane_df,
        left_on="gene_name_prey",
        right_on="gene_name",
        how="left",
        suffixes=("_bait", "_prey")
    )
    del ppi_df["gene_name"]

    return ppi_df


def permute_selection(n_ppi_df, selection_group, target_value, n_permutations=100000, mean=False, permut_fraction=.2):
    category_order = n_ppi_df[selection_group].unique()
    permutation_mat = np.zeros(shape=(n_permutations, len(category_order)))

    sample_size = round((1 - permut_fraction) * n_ppi_df.shape[0])
    for i in range(n_permutations):  # SD on individual proteins not individual PPIs
        ss_n_ppis_df = n_ppi_df.sample(sample_size)
        if mean:
            statistics = ss_n_ppis_df.groupby(selection_group)[target_value].mean()
        else:
            statistics = ss_n_ppis_df.groupby(selection_group)[target_value].sum()
        ordered = statistics.reindex(category_order)
        permutation_mat[i, :] = ordered

    if mean:
        summed_df = n_ppi_df.groupby(selection_group,as_index=False)[target_value].mean()
    else:
        summed_df = n_ppi_df.groupby(selection_group,as_index=False)[target_value].sum()
    summed_df[["ci_2.5", "ci_97.5"]] = np.nanpercentile(permutation_mat,[2.5, 97.5],axis=0).transpose()
    return summed_df

def sum_and_permute(ppi_df, type, dataset, mean):
    gene_column = f"gene_name_{type}"
    membrane_column = f"membrane_association_{type}"
    degree_df = ppi_df.groupby(
        [gene_column, membrane_column],as_index=False).size().dropna().rename({
        membrane_column: "membrane_association"
    },axis=1)
    degree_df = permute_selection(
        degree_df,"membrane_association",
        "size",mean=mean)
    degree_df[["target", "dataset"]] = [type, dataset]
    return degree_df

rule membrane_delta:
    input:
        cvcl_0063_bp=f"work_folder{pn}/data/bioplex/CVCL_0063.csv",
        huri=f"work_folder{pn}/data/huri/intact_huri.csv"
    output:
        membrane_ppis = f"work_folder{pn}/analysis/membrane/HuRI_vs_Bioplex_total_mean.csv"
    run:
        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        huri_df = pd.read_csv(input.huri,sep="\t")

        all_genes = (
                set(bp_df["gene_name_bait"]) |
                set(bp_df["gene_name_prey"]) |
                set(huri_df["gene_name_bait"]) |
                set(huri_df["gene_name_prey"])
        )
        uniprot_id_df = map_membrane_structure(all_genes)
        membrane_df = get_membrane_annotation(uniprot_id_df)
        huri_df = join_membrane_annotation(huri_df,membrane_df)
        bp_df = join_membrane_annotation(bp_df,membrane_df)

        # Mean per dataset
        huri_bait_degree = sum_and_permute(huri_df, "bait", "HuRi", True)
        bp_bait_degree = sum_and_permute(bp_df, "bait", "Bioplex", True)
        huri_prey_degree = sum_and_permute(huri_df, "prey", "HuRi", True)
        bp_prey_degree = sum_and_permute(bp_df, "prey", "Bioplex", True)


        # Shared bait mean
        bait_intersect = set(bp_df["gene_name_bait"]) & set(huri_df["gene_name_bait"])
        shared_bp_df = bp_df[bp_df["gene_name_bait"].isin(bait_intersect)]
        shared_huri_df = huri_df[huri_df["gene_name_bait"].isin(bait_intersect)]


        shared_huri_bait_degree = sum_and_permute(shared_huri_df, "bait", "HuRi", True)
        shared_bp_bait_degree = sum_and_permute(shared_bp_df, "bait", "Bioplex", True)
        shared_huri_prey_degree = sum_and_permute(shared_huri_df, "prey", "HuRi", True)
        shared_bp_prey_degree = sum_and_permute(shared_bp_df, "prey", "Bioplex", True)

        mean_ppis_shared = pd.concat(
            [shared_huri_bait_degree, shared_bp_bait_degree, shared_huri_prey_degree, shared_bp_prey_degree]
        )
        mean_ppis_shared["selection"] = "Shared baits"
        mean_ppis_total = pd.concat(
            [bp_bait_degree, huri_bait_degree, bp_prey_degree, huri_prey_degree]
        )
        mean_ppis_total["selection"] = "All baits"

        pd.concat([mean_ppis_total, mean_ppis_shared]).to_csv(output.membrane_ppis, sep="\t", index=False)


