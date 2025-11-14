import pandas as pd
from mean_distance_support import get_cumulative_sum

rule get_hyrdophobicity_delta:
    """
    Join TASA and THSA from NetSurfP2.0 to POD data
    """
    input:
        uniprot_gene=f"work_folder{pn}/intact/uniprot_to_gene_name.csv",
        rhsa_pdb="data/hydrophobicity/NSP2_complete.tab",
        pod_data=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        hydro_annotated=f"work_folder{pn}/analysis/hydrophobicity/pairs_{{data}}_netsurfp2.csv"
    run:
        ## NetSurfP2.0
        uniprot_2_gene = pd.read_csv(input.uniprot_gene,sep="\t")
        rsha_df = pd.read_csv(input.rhsa_pdb,sep="\t")
        rsha_df = rsha_df.rename({"id": "uniprot_id"},axis=1)
        full = rsha_df.merge(uniprot_2_gene,on="uniprot_id")
        full = full.loc[full.groupby('gene_name')[
            'length'].idxmax()].reset_index(drop=True)  #remove duplicated keeping those with longest length, why not?
        full["gn"] = full["gene_name"]
        hydro_cols = ["thsa_netsurfp2", "tasa_netsurfp2", "rhsa_netsurfp2"]
        full = full[["gn"] + hydro_cols]

        pod_df = pd.read_csv(input.pod_data,sep="\t")
        input_columns = pod_df.columns.values
        pod_df = pod_df.merge(
            full,left_on="gene_name_bait",right_on="gn"
        ).merge(
            full,left_on="gene_name_prey",right_on="gn",suffixes=["_bait", "_prey"]
        )
        pod_df = pod_df.drop(["gn_bait", "gn_prey"],axis=1)
        for c in hydro_cols:
            pod_df[f"{c}_delta"] = abs(pod_df[f"{c}_bait"] - pod_df[f"{c}_prey"])

        all_columns = pod_df.columns.values
        selected_columns = [c for c in all_columns if c not in input_columns]
        pod_df[["pair_id"] + selected_columns].to_csv(output.hydro_annotated,sep="\t",index=False)


rule get_hydro_accumulation:
    """
    Get sliding average of hydrophobicity given POD 
    """
    input:
        pod_data=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
        hydro_pod_data=f"work_folder{pn}/analysis/hydrophobicity/POD_{{data}}_netsurfp2.csv"
    output:
        hydro_lesser=f"work_folder{pn}/analysis/hydrophobicity/cumulative/POD_{{data}}_netsurfp2_lesser.csv",
        hydro_greater=f"work_folder{pn}/analysis/hydrophobicity/cumulative/POD_{{data}}_netsurfp2_greater.csv"
    run:
        measurement_columns = [
            "thsa_netsurfp2_delta",
            "tasa_netsurfp2_delta",
            "rhsa_netsurfp2_delta"
        ]

        hydro_data = pd.read_csv(
            input.hydro_pod_data,sep="\t"
        )
        pod_df = pd.read_csv(
            input.pod_data,sep="\t"
        )

        hydro_pod_df = pod_df.merge(hydro_data, on="pair_id")

        get_cumulative_sum(
            hydro_pod_df,
            value_column="lower_bound_pod",
            cumulative_columns=measurement_columns
        ).to_csv(
            output.hydro_greater,
            sep="\t",index=False
        )

        get_cumulative_sum(
            hydro_pod_df,
            value_column="upper_bound_pod",
            cumulative_columns=measurement_columns,
            greater=False
        ).to_csv(
            output.hydro_lesser,
            sep="\t",index=False
        )
