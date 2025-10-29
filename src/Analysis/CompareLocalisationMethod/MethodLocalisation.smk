import pandas as pd

from ..Annotation.localisation_support import add_localisation
from scipy.stats import wilcoxon, ttest_ind
import seaborn as sns

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
        ""
    run:
        bp_df = pd.read_csv(input.cvcl_0063_bp,sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        huri_df = pd.read_csv(input.huri,sep="\t")
        df_localisation = pd.read_csv(input.localisation_csv,sep="\t")

        bp_df = add_localisation(bp_df,df_localisation)
        huri_df = add_localisation(huri_df,df_localisation)

        bait_intersect = set(bp_df["gene_name_bait"]) & set(huri_df["gene_name_bait"])
        prey_intersect = set(bp_df["gene_name_prey"]) & set(huri_df["gene_name_prey"])
        bp_df = bp_df[bp_df["gene_name_bait"].isin(bait_intersect)]
        huri_df = huri_df[huri_df["gene_name_bait"].isin(bait_intersect)]
        bp_df = bp_df[bp_df["gene_name_prey"].isin(prey_intersect)]
        huri_df = huri_df[huri_df["gene_name_prey"].isin(prey_intersect)]


        # Prey localisation PPI count:
        local_bait_huri = huri_df.groupby(["gene_name_bait", "localisation_bait"],as_index=False).size().rename({
            "size": "huri_ppi"},axis=1)
        local_bait_bp = bp_df.groupby(["gene_name_bait", "localisation_bait"],as_index=False).size().rename({
            "size": "bp_ppi"},axis=1)

        local_bait_ppis = local_bait_bp.merge(local_bait_huri,on=["gene_name_bait", "localisation_bait"],how="outer").fillna(0)
        local_bait_ppis["delta"] = local_bait_ppis["bp_ppi"] - local_bait_ppis["huri_ppi"]
        keep_local = local_bait_ppis.groupby("localisation_bait").size()
        keep_local = keep_local[keep_local >10].index.values

        local_bait_ppis = local_bait_ppis[local_bait_ppis["localisation_bait"].isin(keep_local)]
        sns.kdeplot(data=local_bait_ppis, y="delta", hue="localisation_bait")


        for c_local in local_bait_ppis["localisation_bait"].unique():
            hu_ppi = local_bait_ppis[local_bait_ppis["localisation_bait"] == c_local]["huri_ppi"]
            bp_ppi = local_bait_ppis[local_bait_ppis["localisation_bait"] == c_local]["bp_ppi"]
            s, p = ttest_ind(hu_ppi,bp_ppi)
            print(f"{c_local}, p: {p}")


        # PPI count non overlapping
        bp_df["ppi_id"] = bp_df[["gene_name_bait", "gene_name_prey"]].apply(lambda x: ":".join(x),axis=1)
        huri_df["ppi_id"] = huri_df[["gene_name_bait", "gene_name_prey"]].apply(lambda x: ":".join(x),axis=1)

        shared_ppi = set(bp_df["ppi_id"]) & set(huri_df["ppi_id"])
        shared_prey = set(bp_df["gene_name_prey"]) & set(huri_df["gene_name_prey"])

        bp_df = bp_df[~bp_df["ppi_id"].isin(shared_ppi)]
        huri_df = huri_df[~huri_df["ppi_id"].isin(shared_ppi)]

        localisation_bp = bp_df.groupby(["localisation_bait", "localisation_prey"],as_index=False).size().rename({
            "size": "bp_count"},axis=1)
        localisation_huri = huri_df.groupby(["localisation_bait", "localisation_prey"],as_index=False).size().rename({
            "size": "huri_count"},axis=1)

        localisation_prey = localisation_bp.merge(localisation_huri,on=["localisation_bait", "localisation_prey"],how="outer").fillna(0)
        for c_local in localisation_prey["localisation_bait"].unique():
            hu_ppi = localisation_prey[localisation_prey["localisation_prey"] == c_local]["huri_count"]
            bp_ppi = localisation_prey[localisation_prey["localisation_prey"] == c_local]["bp_count"]
            s, p = wilcoxon(x=hu_ppi-bp_ppi)
            mu = round(np.mean(hu_ppi-bp_ppi),2)
            sd = round(np.std(hu_ppi-bp_ppi),2)
            print(f"mu:{mu}\t sd:{sd}\t p: {p}\t{c_local}")
