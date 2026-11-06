import pandas as pd
import numpy as np


rule get_uniprot_dna_binding:
    output:
        dna_binding = "work_folder/analysis/ProteinPromiscuity/autoactivators/uniprot_tf_activate.tsv"
    run:
        query = "https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true&query=organism_id:9606%20AND%20reviewed:true"
        go_coactive = set(
            pd.read_csv(f"{query}%20AND%20go:0001228&fields=accession", sep="\t", compression="gzip")["Entry"])
        go_activetf = set(
            pd.read_csv(f"{query}%20AND%20go:0003713&fields=accession", sep="\t", compression="gzip")["Entry"])

        df = pd.DataFrame({"uniprot_id": list(go_coactive.union(go_activetf))})
        df.to_csv(output.dna_binding, sep="\t", index=False)


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


rule potential_crapome_proteins:
    params:
        min_tested = 300,
        top = .1
    input:
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/ms_row_wise_detection.tsv",
        crapome_proteins = "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    output:
        crapome_like_proteins = "work_folder/analysis/ProteinPromiscuity/autoactivators/ms_crapome_like_proteins.tsv"
    run:
        df = pd.read_csv(input.row_wise_data, sep="\t")
        crapome_df = pd.read_csv(input.crapome_proteins, sep="\t")
        top_n_proteins = set(crapome_df[crapome_df.mean_sc >= crapome_df.mean_sc.quantile(1 - params.top)]["uniprot_id"])

        bp = df.groupby(["bait", "prey"], as_index=False).agg(detection_rate=("detection", "mean"))
        df_detection_var = bp.groupby("prey", as_index=False).agg(
            detection_var=("detection_rate", "var"),
            detection_mean=("detection_rate", "mean"),
            n_tested=("detection_rate", "size"),
        )
        df_detection_var = df_detection_var[df_detection_var["n_tested"] >= params.min_tested]
        crap_var_dist = df_detection_var[df_detection_var["prey"].isin(top_n_proteins)]["detection_var"]
        crap_mean_dist = df_detection_var[df_detection_var["prey"].isin(top_n_proteins)]["detection_mean"]

        crap_var_quant = crap_var_dist.quantile(0.2)
        crap_mean_quant = crap_mean_dist.quantile(0.2)

        selected_proteins = df_detection_var[
            (df_detection_var["detection_var"] <= crap_var_quant) &
            (df_detection_var["detection_mean"] >= crap_mean_quant)
        ].copy()
        selected_proteins["in_crapome"] = np.where(
            selected_proteins["prey"].isin(crapome_df["uniprot_id"]), "yes", "no")
        selected_proteins.loc[
            selected_proteins["prey"].isin(top_n_proteins), "in_crapome"] = "top"

        selected_proteins.to_csv(output.crapome_like_proteins, sep="\t", index=False)




rule detect_excess_kurtosis:
    params:
        min_pairs = 300
    input:
        ms_adjusted_detectability = "work_folder/analysis/ProteinPromiscuity/detectability/ms_adjusted_detectability.tsv",
        y2h_adjusted_detectability = "work_folder/analysis/ProteinPromiscuity/detectability/y2h_adjusted_detectability.tsv"
    output:
        plot_distribtuions = "work_folder/analysis/ProteinPromiscuity/plot/adjusted_detectability_kurtosis.png"
    threads:
        1
    script:
        "scripts/plot_excess_kurtosis.R"


def what_known(wc):
    if wc.dataset == "ms" and wc.role == "prey":
        return "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    elif wc.dataset == "y2h" and wc.role == "bait":
        return "work_folder/analysis/ProteinPromiscuity/autoactivators/uniprot_tf_activate.tsv"
    else:
        raise ValueError(f"no problematic set for dataset={wc.dataset}, role={wc.role}")


rule p_val_vs_auto:
    input:
        partner_effect = "work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
        selected_proteins = lambda wc: what_known(wc)
    output:
        enrichment="work_folder/analysis/ProteinPromiscuity/autoactivators/known/{dataset}_{role}_added.pdf",
        histogram="work_folder/analysis/ProteinPromiscuity/autoactivators/known/{dataset}_{role}_histogram.pdf"
    script:
        "scripts/check_excess_known_problematic.py"


rule rank_partner_effect:
    input:
        row_wise_data="work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
        known_status = lambda wc: what_known(wc)
    output:
        ranked="work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
    params:
        min_pairs=20,
        n_simulations=1000,
    threads:
        10
    script:
        "scripts/estimate_partner_effect.py"
