from urllib.request import Request

import pandas as pd
import numpy as np
from urllib3.exceptions import RequestError

from ..Annotation.localisation_support import add_localisation
import requests

rule get_bioplex:
    output:
        cvcl_0291_bp="work_folder/data/bioplex/CVCL_0291.csv",
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063.csv"
    shell:
        """
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_293T_DirectedEdges.tsv -O {output.cvcl_0063_bp}
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_HCT116_DirectedEdges.tsv -O {output.cvcl_0291_bp}
        """

rule get_huri:
    """
    must be annotated on gene name as bioplex is reported on gene name 
    """
    input:
        intact="work_folder/formated/bait_prey_publications.csv"
    output:
        huri="work_folder/data/huri/intact_huri.csv"
    run:
        intact = pd.read_csv(input.intact,sep="\t")
        huri_df = intact[intact["pubmed_id"] == 32296183]
        huri_df.to_csv(output.huri,sep="\t",index=False)

rule localisation_delta:
    input:
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063.csv",
        huri="work_folder/data/huri/intact_huri.csv",
        localisation_csv="work_folder/analysis/localisation/gene_to_localisation.csv"
    output:
        localisation_method="work_folder/analysis/localisation/HuRI_vs_Bioplex.csv"
    run:
        n_permutations = 500000

        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        huri_df = pd.read_csv(input.huri,sep="\t")
        df_localisation = pd.read_csv(input.localisation_csv,sep="\t")

        bp_df = add_localisation(bp_df,df_localisation)
        huri_df = add_localisation(huri_df,df_localisation)

        bait_intersect = set(bp_df["gene_name_bait"]) & set(huri_df["gene_name_bait"])
        prey_intersect = set(bp_df["gene_name_prey"]) & set(huri_df["gene_name_prey"])
        shared_baits_bp_df = bp_df[bp_df["gene_name_bait"].isin(bait_intersect)]
        shared_baits_huri_df = huri_df[huri_df["gene_name_bait"].isin(bait_intersect)]


        # Prey localisation PPI count:
        local_bait_huri = shared_baits_huri_df.groupby(["gene_name_bait", "localisation_bait"],as_index=False).size().rename({
            "size": "huri_ppi"},axis=1)
        local_bait_bp = shared_baits_bp_df.groupby(["gene_name_bait", "localisation_bait"],as_index=False).size().rename({
            "size": "bp_ppi"},axis=1)

        local_bait_ppis = local_bait_bp.merge(local_bait_huri,on=["gene_name_bait", "localisation_bait"],how="outer").fillna(0)

        local_bait_ppis["delta"] = local_bait_ppis["huri_ppi"] - local_bait_ppis["bp_ppi"]
        keep_local = local_bait_ppis.groupby("localisation_bait").size()
        keep_local = keep_local[keep_local >10].index.values

        local_bait_ppis = local_bait_ppis[local_bait_ppis["localisation_bait"].isin(keep_local)]

        per_mat = np.zeros(shape=(n_permutations, len(keep_local)))
        sample_size = round((1-0.10)*local_bait_ppis.shape[0])
        for i in range(n_permutations):  # SD on baits not individual ppis
            ss_bait_ppis = local_bait_ppis.sample(sample_size)
            per_mat[i,:] = ss_bait_ppis.groupby("localisation_bait")["delta"].sum().values

        bait_summed_delta = local_bait_ppis.groupby("localisation_bait", as_index=False)["delta"].sum()
        bait_summed_delta["localisation"] = bait_summed_delta["localisation_bait"]
        bait_summed_delta[["ci_2.5", "ci_97.5"]] = np.percentile(per_mat,[2.5, 97.5], axis=0).transpose()
        bait_summed_delta["role"] = "Bait"

        # PPI count non overlapping
        shared_baits_bp_df["ppi_id"] = shared_baits_bp_df[["gene_name_bait", "gene_name_prey"]].apply(lambda x: ":".join(x),axis=1)
        shared_baits_huri_df["ppi_id"] = shared_baits_huri_df[["gene_name_bait", "gene_name_prey"]].apply(lambda x: ":".join(x),axis=1)

        shared_ppi = set(shared_baits_bp_df["ppi_id"]) & set(shared_baits_huri_df["ppi_id"])
        shared_prey = set(shared_baits_bp_df["gene_name_prey"]) & set(shared_baits_huri_df["gene_name_prey"])

        no_overlap_bp_df = shared_baits_bp_df[~shared_baits_bp_df["ppi_id"].isin(shared_ppi)]
        no_overlap_huri_df = shared_baits_huri_df[~shared_baits_huri_df["ppi_id"].isin(shared_ppi)]


        localisation_bp = no_overlap_bp_df.groupby(["gene_name_prey","localisation_prey"],as_index=False).size().rename({
            "size": "bp_ppi"},axis=1)
        localisation_huri = no_overlap_huri_df.groupby(["gene_name_prey","localisation_prey"],as_index=False).size().rename({
            "size": "huri_ppi"},axis=1)
        localisation_prey = localisation_bp.merge(localisation_huri,on=["gene_name_prey","localisation_prey"],how="outer").fillna(0)

        keep_local = localisation_prey.groupby("localisation_prey").size()
        keep_local = keep_local[keep_local >10].index.values
        localisation_prey = localisation_prey[localisation_prey["localisation_prey"].isin(keep_local)]
        localisation_prey["delta"] = localisation_prey["huri_ppi"] - localisation_prey["bp_ppi"]

        per_mat = np.zeros(shape=(n_permutations, len(keep_local)))
        sample_size = round((1-0.10)*localisation_prey.shape[0])
        for i in range(n_permutations):  # SD on baits not individual ppis
            ss_bait_ppis = localisation_prey.sample(sample_size)
            prey_delta = ss_bait_ppis.groupby("localisation_prey")["delta"].sum()
            ordered = prey_delta.reindex(keep_local)
            per_mat[i, :] = ordered

        prey_summed_delta = localisation_prey.groupby("localisation_prey", as_index=False)["delta"].sum()
        prey_summed_delta["localisation"] = prey_summed_delta["localisation_prey"]
        prey_summed_delta[["ci_2.5", "ci_97.5"]] = np.nanpercentile(per_mat,[2.5, 97.5], axis=0).transpose()
        prey_summed_delta["role"] = "Prey"
        columns = ['delta', 'ci_2.5', 'ci_97.5', 'role', 'localisation']
        pd.concat(
            [prey_summed_delta[columns], bait_summed_delta[columns]]
        ).to_csv(output.localisation_method, sep="\t", index=None)


rule membrane_delta:
    input:
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063.csv",
        huri="work_folder/data/huri/intact_huri.csv"
    output:
        localisation_method="work_folder/analysis/localisation/HuRI_vs_Bioplex.csv"
    run:
        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        huri_df = pd.read_csv(input.huri,sep="\t")
        url = "https://rest.uniprot.org/idmapping/run"
        params = {
            "from": "Gene_Name",
            "to": "UniProtKB",
            "ids": ",".join(gene_names),
            "taxonId": "9606"  # Homo sapiens taxonomy ID
        }
        response = requests.post(url,data=params)





def map_membrane_structure(gene_names):
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "Gene_Name",
        "to": "UniProtKB",
        "ids": ",".join(gene_names),
        "taxonId": "9606"  # Homo sapiens taxonomy ID
    }
    response = requests.post(url,data=params)
    if not response.ok:
        raise RequestError("Couldn't submit to uniprot")

    jobid = response.json()["jobId"]

    results_response = requests.get(
        f"https://rest.uniprot.org/idmapping/uniprotkb/results/{jobid}")

