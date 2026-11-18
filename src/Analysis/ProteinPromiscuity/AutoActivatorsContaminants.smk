import pandas as pd
import numpy as np


def dna_binding_info(accessions):
    #TODO check
    out = {}
    for acc in accessions:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
        r.raise_for_status()
        e = r.json()
        feats = [f for f in e.get("features", []) if f["type"] == "DNA binding"]
        kws = {k["id"] for k in e.get("keywords", [])}
        go = {x["id"] for x in e.get("uniProtKBCrossReferences", [])
              if x["database"] == "GO"}
        out[acc] = {
            "dna_bind_feature": [(f["location"]["start"]["value"],
                                  f["location"]["end"]["value"]) for f in feats],
            "keyword_dna_binding": "KW-0238" in kws,
            "go_dna_binding": "GO:0003677" in go,
        }
    return out

print(dna_binding_info(["P04637", "P01308"]))  # TP53 (yes), insulin (no)


rule get_potential_auto_activators:
    params:
        min_tested = 300
    input:
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/y2h_row_wise_detection.tsv",
        auto_activators = "data/autoactivators/y2h.csv"
    output:
        y2h_like_aa = "work_folder/analysis/ProteinPromiscuity/autoactivators/y2h_aa_like_proteins.tsv"
    run:
        df = pd.read_csv(input.row_wise_data, sep="\t")
        bait_proteins = df.bait.unique()

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



rule rank_partner_effect:
    input:
        row_wise_data="work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
        selected_proteins = lambda wc: what_known(wc)
    output:
        ranked="work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
    params:
        min_pairs=20,
        n_simulations=3000,
    run:
        role = wildcards.role
        partner = "prey" if role == "bait" else "bait"
        df = pd.read_csv(input.row_wise_data, sep="\t")
        pair = df.groupby([role, partner], as_index=False).agg(
            n=("detection", "size"), k=("detection", "sum")
        )
        pair = pair[pair["n"] >= 2].copy()     # n==1 pairs are uninformative by construction
        pair["det"] = pair["k"] > 0            # detected at least once
        pair["full"] = pair["k"] == pair["n"]  # detected on every test
        tot = pair.groupby(role).agg(
            n_pairs=("n", "size"), n_obs=("n", "sum"), n_pos=("k", "sum"),
            n_det=("det", "sum"), n_full=("full", "sum"),
        )
        tot["p"] = tot["n_pos"] / tot["n_obs"]
        pair = pair.join(tot["p"], on=role)
        e = pair["n"] * pair["p"]
        v = pair["n"] * pair["p"] * (1 - pair["p"])
        pair["chi"] = ((pair["k"] - e) ** 2 / v).where(v > 0, 0.0)
        tot["D"] = pair.groupby(role)["chi"].sum() / (tot["n_pairs"] - 1)
        tot["f"] = (tot["n_full"] / tot["n_det"]).where(tot["n_det"] > 0) # Fraction of full replication

        # p = 0 or 1: every pair is uniform by construction, so D carries no information
        tot.loc[(tot["p"] == 0) | (tot["p"] == 1), "D"] = np.nan
        tot = tot[tot["n_pairs"] >= params.min_pairs]

        rng = np.random.default_rng(0)

        def uniformity_p(g):
            n, p, r, B = g["n"].to_numpy(), g["p"].iloc[0], params.rho0, params.n_simulations
            sim = rng.binomial(n, rng.beta(p * (1 - r) / r, (1 - p) * (1 - r) / r, size=(B, n.size)))
            q = sim.sum(axis=1, keepdims=True) / n.sum()  # re-estimate p as the observed stat does
            chi = ((sim - n * q) ** 2 / np.maximum(n * q * (1 - q), 1e-12)).sum(axis=1)
            return (1 + (chi <= g["chi"].sum()).sum()) / (B + 1)

        tot["p_vs_specific"] = (
            pair[pair[role].isin(tot.index[tot["D"].notna()])].groupby(role).apply(uniformity_p)
        )

        detected = tot["n_pos"] > 0
        tot["score"] = np.nan
        tot.loc[detected, "score"] = (-tot.loc[detected, "D"]).rank(
            pct=True
        ) + tot.loc[detected, "p"].rank(pct=True)
        tot.sort_values("score", ascending=False, na_position="last").to_csv(
            output.ranked, sep="\t"
        )
