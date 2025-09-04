import pandas as pd

#TODO: generalise this for all similarity

def get_sliding_avg_enrichment_hydro(df, value_column, greater=True, min_samples=50):
    bins = df[value_column].unique()
    bins.sort()
    values = df[value_column].values
    idx_val = values.argsort()
    if greater:
        idx_val = idx_val[::-1]
        bins = bins[::-1]
    values = values[idx_val]

    thsa = df["thsa_netsurfp2"].to_numpy()[idx_val]
    tasa = df["tasa_netsurfp2"].to_numpy()[idx_val]
    rhsa = df["rhsa_netsurfp2"].to_numpy()[idx_val]

    if not greater:
        bins = -bins
        values = -values

    thsa_sum = 0
    tasa_sum = 0
    rhsa_sum = 0

    previous = 0
    i = 0
    j = 0
    rows = [{}] * len(bins)
    for threshold in bins:
        while (i < len(values) and threshold <= values[i]) or i < min_samples:
            i += 1


        thsa_sum += thsa[previous:i].sum()
        tasa_sum += tasa[previous:i].sum()
        rhsa_sum += rhsa[previous:i].sum()

        rows[j] = {
            "limit": value_column,
            "value": (threshold if greater else -threshold),
            "thsa_avg": thsa_sum / i,
            "tasa_avg": tasa_sum / i,
            "rhsa_avg": rhsa_sum / i,
            "number_of_pairs": i
        }
        if previous != i:
            j += 1
        previous = i

    return pd.DataFrame([r for r in rows if r])

rule get_hyrdophobicity_delta:
    input:
        uniprot_gene="work_folder/intact/uniprot_to_gene_name.csv",
        rhsa_pdb="data/hydrophobicity/NSP2_complete.tab",
        flat_file="work_folder/inferred_search_space/analysis/bias_reduced_ppis/p_estimated_protein_pairs.csv",
        abundance_file="work_folder/analysis/abundance_aware/POD_abundance.csv"
    output:
        abundance_out="work_folder/analysis/hydrophobicity/abunadce_netsurfp2.0.csv",
        flat_out="work_folder/analysis/hydrophobicity/flat_netsurfp2.0.csv"
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

        ## Abundance
        abundance_df = pd.read_csv(input.abundance_file,sep="\t")
        abundance_df = abundance_df.merge(
            full,left_on="gene_name_bait",right_on="gn"
        ).merge(
            full,left_on="gene_name_prey",right_on="gn",suffixes=["_bait", "_prey"]
        )
        abundance_df = abundance_df.drop(["gn_bait", "gn_prey"],axis=1)
        for c in hydro_cols:
            abundance_df[f"{c}_delta"] = abs(abundance_df[f"{c}_bait"] - abundance_df[f"{c}_prey"])

        abundance_df.to_csv(output.abundance_out,sep="\t",index=False)

        ## Flat
        flat_df = pd.read_csv(input.flat_file,sep="\t")
        flat_df = flat_df.merge(
            full,left_on="gene_name_bait",right_on="gn"
        ).merge(
            full,left_on="gene_name_prey",right_on="gn",suffixes=["_bait", "_prey"]
        )
        flat_df = flat_df.drop(["gn_bait", "gn_prey"],axis=1)
        for c in hydro_cols:
            flat_df[f"{c}_delta"] = abs(flat_df[f"{c}_bait"] - flat_df[f"{c}_prey"])

        flat_df.to_csv(output.flat_out,sep="\t",index=False)

rule get_go_accumulation:
    input:
        abundance_in="work_folder/analysis/hydrophobicity/abunadce_netsurfp2.0.csv",
        flat_in="work_folder/analysis/hydrophobicity/flat_netsurfp2.0.csv"
    output:
        flat_netsurfp_greater="work_folder/analysis/hydrophobicity/flat_netsurfp_greater.csv",
        abundance_netsurfp_greater="work_folder/analysis/hydrophobicity/abundance_netsurfp_greater.csv",
        flat_netsurfp_lesser="work_folder/analysis/hydrophobicity/flat_netsurfp_lesser.csv",
        abundance_netsurfp_lesser="work_folder/analysis/hydrophobicity/abundance_netsurfp_lesser.csv"
    run:
        abundance_df = pd.read_csv(
            input.abundance_in, sep="\t"
        )
        get_sliding_avg_enrichment_hydro(
            abundance_df,
            "lower_bound_pod").to_csv(
            output.abundance_netsurfp_greater,
            sep="\t",index=False
        )
        get_sliding_avg_enrichment_hydro(
            abundance_df,
            "upper_bound_pod",
            greater=False).to_csv(
            output.abundance_netsurfp_lesser,
            sep="\t", index=False
        )

        flat_df = pd.read_csv(
            input.flat_in, sep="\t"
        )
        get_sliding_avg_enrichment_hydro(
            flat_df,
            "p_lower_ci").to_csv(
            output.flat_netsurfp_greater,
            sep="\t",index=False
        )
        get_sliding_avg_enrichment_hydro(
            flat_df,
            "p_upper_ci",
            greater=False).to_csv(
            output.flat_netsurfp_lesser,
            sep="\t", index=False
        )
