import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from py_qvalue import qvalue

prot_column = "protein"
ROLE_COLORS = {"bait": "tab:blue", "prey": "tab:orange"}


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
    features = df["feature"].unique()
    x = np.arange(len(features))
    width = 0.2
    fig, ax = plt.subplots(figsize=(7, 4))
    # per feature: [bait sig, bait other, prey sig, prey other]
    for i, (role, color) in enumerate(ROLE_COLORS.items()):
        role_df = df[df["role"] == role].set_index("feature").loc[features]
        offset = (2 * i - 1.5) * width
        ax.bar(x + offset, role_df["sig_r"], width, color=color, edgecolor="black")
        ax.bar(x + offset + width, role_df["other_r"], width, color=color, edgecolor="black", hatch="//")

    ax.set_xticks(x)
    ax.set_xticklabels(features)
    ax.set_xlabel("feature")
    ax.set_ylabel("value")

    ax.legend(
        handles=[Patch(facecolor=color, edgecolor="black", label=role) for role, color in ROLE_COLORS.items()]
        + [
            Patch(facecolor="white", edgecolor="black", label="sig_r"),
            Patch(facecolor="white", edgecolor="black", hatch="//", label="other_r"),
        ]
    )

    plt.tight_layout()
    plt.savefig(snakemake.output.plot_ratio, format="svg", bbox_inches="tight")


def plot_or(df):
    order = df.groupby("feature")["odds_ratio"].mean().sort_values().index
    y = pd.Series(np.arange(len(order)), index=order)

    fig, ax = plt.subplots(
        figsize=(6, 0.5 * len(order) + 1.5))
    for i, (role, color) in enumerate(ROLE_COLORS.items()):
        role_df = df[df["role"] == role]
        ax.errorbar(
            role_df["odds_ratio"],
            y[role_df["feature"]].to_numpy() + (i - 0.5) * 0.3,
            xerr=[role_df["odds_ratio"] - role_df["odds_low"], role_df["odds_high"] - role_df["odds_ratio"]],
            fmt="o",
            capsize=4,
            color=color,
            label=role,
        )
    ax.axvline(1, color="grey", ls="--")
    ax.set(
        xscale="log", yticks=y.to_numpy(), yticklabels=order, xlabel="Odds ratio (95% CI)"
    )
    ax.legend()

    plt.tight_layout()
    plt.savefig(snakemake.output.plot_or, bbox_inches="tight")
    

def read_role_effects(path):
    df = pd.read_csv(path, sep="\t", index_col=0).dropna(subset=["p_vs_specific"])
    role = df.index.name  # estimate_partner_effect.py writes the role as the index name
    df = df.rename_axis(prot_column).reset_index()
    df["role"] = role
    # q-values within role, as each role is its own set of tests
    df["quvalue"] = qvalue(df["p_vs_specific"].to_numpy(), lfdr_out=False)["qvalues"]
    df["sig"] = df["quvalue"] < 0.05
    return df


def main():
    df_protein_effects = pd.concat(
        [read_role_effects(path) for path in snakemake.input.ranked], ignore_index=True
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

    feature_res = []
    for role, role_df in df_protein_effects.groupby("role"):
        for feature in features_to_test:
            annotated = role_df[role_df[feature].notna()]
            sig = annotated["sig"]
            sig_r = annotated.loc[sig, feature].mean()
            other_r = annotated.loc[~sig, feature].mean()
            cross_tab = pd.crosstab(annotated["sig"], annotated[feature])
            c_or, p = fisher_exact(cross_tab)
            res = odds_ratio(cross_tab, kind="sample")
            ci = res.confidence_interval(confidence_level=0.95)
            feature_res.append(
                [
                    role,
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
            "role",
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
