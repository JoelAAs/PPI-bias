import pandas as pd
import numpy as np


rule get_top_crap:
    input:
        crapome = "data/autoactivators/CrapOme6926.csv"
    output:
        crapome_proteins = "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    run:
        df = pd.read_csv(input.crapome, sep="\t")
        df.set_index("Uniprot ID", inplace=True)
        drop = ["geneSymbol", "RefSeq Protein ID", "Ensemble Protein ID", "Entrez Gene ID"]
        df.drop(columns=drop, inplace=True)
        
        df["mean_sc"] = df.mean(axis=1)
        
        df["uniprot_id"] = df.index
        df = df[df.mean_sc > 1]
        df[["uniprot_id", "mean_sc"]].to_csv(output.crapome_proteins, sep="\t", index=False)


rule rank_partner_effect:
    input:
        row_wise_data="work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
    output:
        ranked="work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
    params:
        min_pairs=20,
        n_simulations=1000,
    threads:
        10
    script:
        "scripts/estimate_partner_effect.py"


rule get_abundance_paxdb:
    input:
        idmapping = "work_folder/data/uniprot/HUMAN_9606_idmapping.dat.gz"
    output:
        pax_db_data="work_folder/data/paxdb_human.tsv"
    script:
        "scripts/get_format_paxdb.py"


rule get_localisation_file:
    input:
        localisation=config["localisation_file"],
        idmapping="work_folder/data/uniprot/HUMAN_9606_idmapping.dat.gz"
    output:
        with_uniprot="work_folder/data/localisation/uniprot_localisation.csv"
    run:
        localisation_df = pd.read_csv(input.localisation, sep="\t")
        localisation_df = localisation_df[localisation_df["Reliability"].isin(["Enhanced", "Supported", "Approved"])]

        idmapping = pd.read_csv(
            input.idmapping, sep="\t", header=None,
            names=["uniprot_id", "type", "Gene"]
        )
        idmapping = idmapping[idmapping["type"] == "Ensembl"]
        idmapping["Gene"] = idmapping["Gene"].str.replace(r"\.\d+$", "", regex=True)

        # accessions mapped to genes with different main locations are dropped
        loc_df = localisation_df.merge(idmapping[["uniprot_id", "Gene"]], on="Gene")[["uniprot_id", "Main location"]] \
            .drop_duplicates() \
            .drop_duplicates("uniprot_id", keep=False)
        loc_df.to_csv(output.with_uniprot, sep="\t", index=False)


rule group_uniform:
    input:
        ranked="work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
        pax_db_data="work_folder/data/paxdb_human.tsv",
        localisation="work_folder/data/localisation/uniprot_localisation.csv",
        crapome_proteins="work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    output:
        annotated="work_folder/analysis/ProteinPromiscuity/autoactivators/groups/{dataset}_{role}_partner_effect.tsv",
        tests="work_folder/analysis/ProteinPromiscuity/autoactivators/groups/{dataset}_{role}_feature_tests.tsv",
        plot_or="work_folder/analysis/ProteinPromiscuity/autoactivators/groups/{dataset}_{role}_feature_log_or.svg",
        plot_ratio="work_folder/analysis/ProteinPromiscuity/autoactivators/groups/{dataset}_{role}_feature_ratio.svg"
    script:
        "scripts/significant_grouping.py"