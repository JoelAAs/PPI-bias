import pandas as pd
from scipy.stats import beta


def get_degree(df, bait=True):
    pseudo_n = 1
    df["ratio"] = df["n_observed"] / df["n_tested"]
    mean_p = df["ratio"].mean()  # probability of an interaction given a random protein pair
    prior_alpha = pseudo_n * mean_p
    prior_beta = pseudo_n - prior_alpha
    p_prior = prior_alpha / (prior_alpha + prior_beta)
    low_prior = beta.ppf(0.025,prior_alpha,prior_beta)
    high_prior = beta.ppf(0.975,prior_alpha,prior_beta)

    suffix = ["bait", "prey"] if bait else ["prey", "bait"]
    id_col = f"gene_name_{suffix[0]}"
    other_col = f"gene_name_{suffix[1]}"

    n_unique = df[other_col].nunique()
    prob_cols = ["p", "lower_bound_pod", "upper_bound_pod"]
    df_degree = df.groupby(id_col,as_index=False)[prob_cols].sum()
    df_tests = df.groupby(id_col,as_index=False).size().rename({id_col: "gene_name"},axis=1)
    df_tests[f"num_untested_{suffix[0]}"] = n_unique - df_tests["size"]
    del df_tests["size"]
    df_degree = df_degree.rename({
        f"gene_name_{suffix[0]}": "gene_name",
        "p": f"obs_mean_{suffix[0]}_degree",
        "lower_bound_pod": f"obs_lower_{suffix[0]}_degree",
        "upper_bound_pod": f"obs_upper_{suffix[0]}_degree"
    },axis=1)
    df_degree = df_degree.merge(df_tests,on="gene_name")

    df_degree[f"mean_{suffix[0]}_degree"] = df_degree[f"obs_mean_{suffix[0]}_degree"] + df_degree[
        f"num_untested_{suffix[0]}"] * p_prior
    df_degree[f"lower_{suffix[0]}_degree"] = df_degree[f"obs_lower_{suffix[0]}_degree"] + df_degree[
        f"num_untested_{suffix[0]}"] * low_prior
    df_degree[f"upper_{suffix[0]}_degree"] = df_degree[f"obs_upper_{suffix[0]}_degree"] + df_degree[
        f"num_untested_{suffix[0]}"] * high_prior

    return df_degree, [p_prior, low_prior, high_prior, n_unique]


def fill_na(df, params_bait, params_prey):
    suffix = ["bait", "prey"]
    df[f"mean_{suffix[0]}_degree"] = df[f"mean_{suffix[0]}_degree"].fillna(params_bait[-1] * params_bait[0])
    df[f"lower_{suffix[0]}_degree"] = df[f"lower_{suffix[0]}_degree"].fillna(params_bait[-1] * params_bait[1])
    df[f"upper_{suffix[0]}_degree"] = df[f"upper_{suffix[0]}_degree"].fillna(params_bait[-1] * params_bait[2])
    df[f"mean_{suffix[1]}_degree"] = df[f"mean_{suffix[1]}_degree"].fillna(params_prey[-1] * params_prey[0])
    df[f"lower_{suffix[1]}_degree"] = df[f"lower_{suffix[1]}_degree"].fillna(params_prey[-1] * params_prey[1])
    df[f"upper_{suffix[1]}_degree"] = df[f"upper_{suffix[1]}_degree"].fillna(params_prey[-1] * params_prey[2])
    df[f"num_untested_{suffix[0]}"] = df[f"num_untested_{suffix[0]}"].fillna(params_bait[-1])
    df[f"num_untested_{suffix[1]}"] = df[f"num_untested_{suffix[1]}"].fillna(params_prey[-1])

    return df


def threshold_degree(df, t, mode="interaction"):
    modes = ["interaction", "non_interaction", "all"]
    if mode not in modes:
        raise ValueError(f"Unknown mode: {mode}")
    if mode == "interaction":
        df_t = df[df["lower_bound_pod"] > t]
    elif mode == "non_interaction":
        df_t = df[(df["n_observed"] == 0) & (df["n_tested"] >= t)]
    else:
        df_t = df[df["n_observed"] != 0]

    df_bait_degree = df_t.groupby("gene_name_bait",as_index=False).size()
    df_bait_degree = df_bait_degree.rename({
        "gene_name_bait": "gene_name",
        "size": "degree_bait"
    },axis=1)
    df_prey_degree = df_t.groupby("gene_name_prey",as_index=False).size()
    df_prey_degree = df_prey_degree.rename({
        "gene_name_prey": "gene_name",
        "size": "degree_prey"
    },axis=1)
    t_degree = df_bait_degree.merge(df_prey_degree,on="gene_name",how="outer").fillna(0)

    return t_degree


rule get_degree_dist_hippie:
    """
    Get degree list from HIPPIE-current
    """
    input:
        hippie="data/HIPPIE-current.mitab.txt"
    output:
        degree=f"work_folder{pn}/degree/full_hippie.csv"
    run:
        df_hippie = pd.read_csv(input.hippie,sep="\t")
        gene_cols = ["Gene Name Interactor A", "Gene Name Interactor B"]
        df_hippie = df_hippie.dropna(subset=gene_cols).loc[
            (df_hippie["Taxid Interactor A"] == "taxid:9606(Homo sapiens)") &
            (df_hippie["Taxid Interactor B"] == "taxid:9606(Homo sapiens)")
            ]

        df_hippie["interaction_id"] = df_hippie[gene_cols].apply(lambda x: "-".join(sorted(x)),axis=1)
        df_hippie_long = pd.melt(df_hippie,id_vars=["interaction_id"],value_vars=gene_cols,value_name="gene_name")
        hippie_degree = df_hippie_long.groupby("gene_name",as_index=False)["interaction_id"].nunique().rename(
            {"interaction_id": "degree"},axis=1
        )
        hippie_degree.to_csv(output.degree,sep="\t",index=False)

rule get_degree_values:
    """
    Get degree distribution for flat POD
    """
    params:
        hci_limits=config["hci_limits"],
        hcni_tested=config["hcni_tested"]
    input:
        pod_file=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv"
    output:
        summed_probability=f"work_folder{pn}/degree/{{data}}_summed.csv",
        naive_degree = f"work_folder{pn}/degree/{{data}}_naive.csv",
        hci_threshold=expand(
            "work_folder{pn}/degree/{{data}}_HCI_{hci_limit}.csv",
            hci_limit=config["hci_limits"],pn=pn),
        hcni_tests=expand(
            "work_folder{pn}/degree/{{data}}_HCNI_{hcni_tested}.csv",
            hcni_tested=config["hcni_tested"],pn=pn)
    run:
        df = pd.read_csv(input.pod_file,sep="\t")

        bait_degree, params_bait = get_degree(df)
        prey_degree, params_prey = get_degree(df,False)
        full_degree = bait_degree.merge(prey_degree,on="gene_name",how="outer")
        full_degree = fill_na(full_degree,params_bait,params_prey)
        full_degree["degree_bait"] = full_degree["lower_bait_degree"]
        full_degree["degree_prey"] = full_degree["lower_prey_degree"]
        full_degree.to_csv(output.summed_probability,sep="\t",index=False)
        for hci_filename, hci_limit in zip(output.hci_threshold,params.hci_limits):
            df_hci = threshold_degree(df,hci_limit)
            df_hci.to_csv(hci_filename,sep="\t",index=False)

        for hcni_filename, hcni_limit in zip(output.hcni_tests,params.hcni_tested):
            df_hcni = threshold_degree(df,hcni_limit,mode="non_interaction")
            df_hcni.to_csv(hcni_filename,sep="\t",index=False)

        threshold_degree(df,0,mode="all").to_csv(output.naive_degree,sep="\t",index=False)

rule get_intact_bait_prey_degree:
    input:
        intact_bp = "work_folder/per_gene/formated/bait_prey_publications.csv"
    output:
        intact_degree = f"work_folder{pn}/degree/intact.csv"
    run:
        df_intact = pd.read_csv(input.intact_bp, sep="\t")
        df_bait_degree = df_intact.groupby("gene_name_bait",as_index=False).size()
        df_bait_degree = df_bait_degree.rename({
            "gene_name_bait": "gene_name",
            "size": "degree_bait"
        },axis=1)
        df_prey_degree = df_intact.groupby("gene_name_prey",as_index=False).size()
        df_prey_degree = df_prey_degree.rename({
            "gene_name_prey": "gene_name",
            "size": "degree_prey"
        },axis=1)
        t_degree = df_bait_degree.merge(df_prey_degree,on="gene_name",how="outer").fillna(0)
        t_degree.to_csv(output.intact_degree,sep="\t",index=False)
