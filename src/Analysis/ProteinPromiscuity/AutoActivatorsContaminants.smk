import pandas as pd
import numpy as np

rule get_potential_auto_activators:
    params:
        min_tested = 300
    input:
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/y2h_row_wise_detection.tsv",
        known_auto_activators = "data/autoactivators/y2h.csv"
    output:
        y2h_like_aa = "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_aa_like_proteins.tsv"
    run:
        df = pd.read_csv(input.row_wise_data, sep="\t")
        y2h_df = pd.read_csv(input.known_auto_activators, sep="\t")
        known_proteins = set(y2h_df["UniProt"])

        bp = df.groupby(["bait", "prey"], as_index=False).agg(detection_rate=("detection", "mean"))
        df_detection_var = bp.groupby("bait", as_index=False).agg(
            detection_var=("detection_rate", "var"),
            detection_mean=("detection_rate", "mean"),
            n_tested=("detection_rate", "size"),
        )
        df_detection_var = df_detection_var[df_detection_var["n_tested"] >= params.min_tested]
        aa_var_dist = df_detection_var[df_detection_var["bait"].isin(known_proteins)]["detection_var"]
        aa_mean_dist = df_detection_var[df_detection_var["bait"].isin(known_proteins)]["detection_mean"]

        aa_var_quant = aa_var_dist.quantile(0.2)
        aa_mean_quant = aa_mean_dist.quantile(0.2)

        selected_proteins = df_detection_var[
            (df_detection_var["detection_var"] <= aa_var_quant) &
            (df_detection_var["detection_mean"] >= aa_mean_quant)
        ].copy()
        selected_proteins["known_auto_activator"] = np.where(
            selected_proteins["bait"].isin(known_proteins), "yes", "no")
        selected_proteins["uniprot_id"] = selected_proteins["bait"]

        selected_proteins.to_csv(output.y2h_like_aa, sep="\t", index=False)


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


# def what_crap(wc):
#     if wc.dataset == "ms":
#         return "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
#     elif wc.dataset == "y2h":
#         return "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_aa_like_proteins.tsv"
#     else:
#         raise Error("wrong wildcard", wc.dataset)

rule model_select_genes:
    input:        
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
        selected_proteins = lambda wc: what_crap(wc)
    output:
        crap_model = "work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_model.jls",
    log:
        "logs/analysis/ProteinPromiscuity/{dataset}_{role}_fit.log",
    conda:
        "julia"
    threads:
        15
    shell:
        """
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} julia -t 1 src/Analysis/ProteinPromiscuity/scripts/fit_selected_model.jl \
        {input.row_wise_data} {input.selected_proteins} {wildcards.role} {output.crap_model} \
        > {log} 2>&1
        """



rule plot_partner_effects:
    params:
        crapome_top = .1
    input:
        ranked_y2h = "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_bait_partner_effect.tsv",
        ranked_ms = "work_folder/analysis/ProteinPromiscuity/autoactivators/ms_prey_partner_effect.tsv",
        known_auto_activators = "data/autoactivators/y2h.csv",
        crapome_proteins = "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    output:
        partner_effect_plot = "work_folder/analysis/ProteinPromiscuity/plot/partner_effect.png",
        validation_plot = "work_folder/analysis/ProteinPromiscuity/plot/partner_effect_validation.png"
    threads:
        1
    script:
        "scripts/plot_partner_effects.R"



rule plot_bait_prey_detectability:
    input:
        ms_bait = "work_folder/analysis/ProteinPromiscuity/autoactivators/ms_bait_partner_effect.tsv",
        ms_prey = "work_folder/analysis/ProteinPromiscuity/autoactivators/ms_prey_partner_effect.tsv",
        y2h_bait = "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_bait_partner_effect.tsv",
        y2h_prey = "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_prey_partner_effect.tsv"
    output:
        detectability = "work_folder/analysis/ProteinPromiscuity/plot/naive_detectability.png",
        ms_adjusted_detectability = "work_folder/analysis/ProteinPromiscuity/detectability/ms_adjusted_detectability.tsv",
        y2h_adjusted_detectability = "work_folder/analysis/ProteinPromiscuity/detectability/y2h_adjusted_detectability.tsv"
    threads:
        1
    script:
        "scripts/plot_bait_prey_detectability.R"



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
    if wc.dataset == "ms":
        return "work_folder/analysis/ProteinPromiscuity/autoactivators/crapome_proteins.tsv"
    elif wc.dataset == "y2h":
        return "data/autoactivators/y2h.csv"
    else:
        raise Error("wrong wildcard", wc.dataset)

rule p_val_vs_auto:
    input:
        partner_effect = "work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
        selected_proteins = lambda wc: what_known(wc)
    output:
        enrichment="work_folder/analysis/ProteinPromiscuity/autoactivators/known/{dataset}_{role}_added.pdf",
        histogram="work_folder/analysis/ProteinPromiscuity/autoactivators/known/{dataset}_{role}_histogram.pdf"
    script:
        "scripts/check_excess_known_problematic.py"

