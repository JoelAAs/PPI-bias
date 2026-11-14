import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from py_qvalue import qvalue

prot_column = snakemake.wildcards.role


# MS
def annotate_abundance(df, high_fraction=0.10):
    paxdb_df = pd.read_csv(
        snakemake.input.pax_db_data, sep="\t", usecols=["uniprot_id", "abundance"]
    ).rename({"uniprot_id": prot_column}, axis=1)
    top_quant_abundance = paxdb_df["abundance"].quantile(1 - high_fraction)
    paxdb_df["abundance"] = paxdb_df["abundance"] >= top_quant_abundance

    df = df.merge(paxdb_df, on=prot_column, how="left")
    return df


def annotate_crapome_status(df):
    df_crapome = pd.read_csv(
        snakemake.input.crapome_proteins,
        sep="\t"
    )
    df["crapome"] = df[prot_column].isin(df_crapome["uniprot_id"])
    return df


def annotate_folding(df):
    query = "https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true&query=organism_id:9606"
    go_chap = set(
        pd.read_csv(
            f"{query}%20AND%20go:0044183&fields=accession", sep="\t", compression="gzip"
        )["Entry"]
    )
    go_hold = set(
        pd.read_csv(
            f"{query}%20AND%20go:0140309&fields=accession", sep="\t", compression="gzip"
        )["Entry"]
    )

    folding_proteins = go_chap | go_hold
    df["folding_proteins"] = df[prot_column].isin(folding_proteins)
    return df


def annotate_nucleotide_binding(df):
    query = "https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true&query=organism_id:9606"
    nulc_binding = set(
        pd.read_csv(
            f"{query}%20AND%20GO:0003676&fields=accession", sep="\t", compression="gzip"
        )["Entry"]
    )

    df["nulc_binding"] = df[prot_column].isin(nulc_binding)
    return df


def annotate_common_contaminant(df):
    # https://doi.org/10.1021/acs.jproteome.2c00145
    contaminants_df = pd.read_excel(
        "https://raw.githubusercontent.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/main/Universal%20protein%20contaminant%20FASTA/Contaminant%20protein%20information%20in%20the%20FASTA%20library%20and%20the%20potential%20source%20of%20contaminations.xlsx",
        skiprows=1,
    )
    df["contaminant"] = df[prot_column].isin(contaminants_df["Uniprot ID"])

    return df


def annotate_ms(df):
    df = annotate_abundance(df)
    df = annotate_folding(df)
    df = annotate_nucleotide_binding(df)
    df = annotate_common_contaminant(df)
    df = annotate_crapome_status(df)
    return df

# y2h
def annotate_nucleoplasm_localisation(df):
    loc_df = pd.read_csv(snakemake.input.localisation, sep="\t").rename(
        {"uniprot_id": prot_column}, axis=1
    )
    loc_df["nucleoplasm"] = loc_df["Main location"].str.contains("Nucleoplasm")
    df = df.merge(loc_df[[prot_column, "nucleoplasm"]], on=prot_column, how="left")
    return df

def annotate_coactive_tf(df):
    query = "https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true&query=organism_id:9606"
    go_coactive = set(
        pd.read_csv(f"{query}%20AND%20go:0001228&fields=accession", sep="\t", compression="gzip")["Entry"])
    go_activetf = set(
        pd.read_csv(f"{query}%20AND%20go:0003713&fields=accession", sep="\t", compression="gzip")["Entry"])
    tfs = go_coactive | go_activetf

    df["coac_tfs"] = df[prot_column].isin(tfs)

    return df

def annotate_y2h(df):
    df = annotate_nucleoplasm_localisation(df)
    df = annotate_folding(df)
    df = annotate_coactive_tf(df)
    return df


def plot_ratio(df):
    x = np.arange(len(df))
    width = 0.4
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x - width / 2, df["sig_r"], width, color="white", edgecolor="black")
    ax.bar(x + width / 2, df["other_r"], width, color="white", edgecolor="black", hatch="//"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(df["feature"])
    ax.set_xlabel("feature")
    ax.set_ylabel("value")

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor="white", edgecolor="black", label="sig_r"),
            Patch(facecolor="white", edgecolor="black", hatch="//", label="other_r"),
        ]
    )

    plt.tight_layout()
    plt.savefig(snakemake.output.plot_ratio, format="svg", bbox_inches="tight")


def plot_or(df):
    df = df.sort_values("odds_ratio")
    y = np.arange(len(df))

    fig, ax = plt.subplots(
        figsize=(6, 0.5 * len(df) + 1.5))
    ax.errorbar(
        df["odds_ratio"],
        y,
        xerr=[df["odds_ratio"] - df["odds_low"], df["odds_high"] - df["odds_ratio"]],
        fmt="o",
        capsize=4,
    )
    ax.axvline(1, color="grey", ls="--")
    ax.set(
        xscale="log", yticks=y, yticklabels=df["feature"], xlabel="Odds ratio (95% CI)"
    )

    plt.tight_layout()
    plt.savefig(snakemake.output.plot_or, bbox_inches="tight")
    

def main():
    df_protein_effects = pd.read_csv(snakemake.input.ranked, sep="\t").dropna(
        subset=["p_vs_specific"]
    )

    dataset = snakemake.wildcards.dataset

    if dataset == "ms":
        df_protein_effects = annotate_ms(df_protein_effects)
    elif dataset == "y2h":
        df_protein_effects = annotate_y2h(df_protein_effects)
    else:
        raise ValueError(f"Unkown dataset: {dataset}")

    features_to_test = {
        "ms": ["abundance", "folding_proteins", "nulc_binding", "contaminant", "crapome"],
        "y2h": ["nucleoplasm", "folding_proteins", "coac_tfs"],
    }[dataset]

    df_protein_effects["quvalue"] = qvalue(
        df_protein_effects["p_vs_specific"].to_numpy(), lfdr_out=False
    )["qvalues"]
    df_protein_effects["sig"] = df_protein_effects["quvalue"] < 0.05

    feature_res = []
    for feature in features_to_test:
        annotated = df_protein_effects[df_protein_effects[feature].notna()]
        sig = annotated["sig"]
        sig_r = annotated.loc[sig, feature].mean()
        other_r = annotated.loc[~sig, feature].mean()
        cross_tab = pd.crosstab(annotated["sig"], annotated[feature])
        c_or, p = fisher_exact(cross_tab)
        res = odds_ratio(cross_tab, kind="sample")
        ci = res.confidence_interval(confidence_level=0.95)
        feature_res.append(
            [
                feature,
                annotated["sig"].sum(),
                (~annotated["sig"]).sum(),
                sig_r,
                other_r,
                c_or,
                ci.low,
                ci.high,
                p,
            ]
        )

    feature_res = pd.DataFrame(
        feature_res,
        columns=[
            "feature",
            "n_sig",
            "n_other",
            "sig_r",
            "other_r",
            "odds_ratio",
            "odds_low",
            "odds_high",
            "p",
        ],
    )

    plot_or(feature_res)
    plot_ratio(feature_res)

    df_protein_effects.to_csv(snakemake.output.annotated, sep="\t", index=False)
    feature_res.to_csv(snakemake.output.tests, sep="\t", index=False)

main()
