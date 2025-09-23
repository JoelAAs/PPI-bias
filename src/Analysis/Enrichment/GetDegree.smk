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
    prob_cols = ["p", "lower_bound_pod", "lower_bound_pod"]
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


def threshold_degree(df, t):
    df_t = df[df["lower_bound_pod"] > t]
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
    input:
        hippie="data/HIPPIE-current.mitab.txt"
    output:
            degree="work_folder/degree/full_hippie.csv"
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
            {"interaction_id": "degree"}, axis=1
        )
        hippie_degree.to_csv(output.degree,sep="\t",index=False)

rule flat_degree_dist:
    input:
        flat_probability="work_folder/analysis/POD/POD_flat.csv"
    output:
        summed_probability="work_folder/degree/flat_summed.csv",
        threshold_1="work_folder/degree/flat_min.1.csv",
        threshold_2="work_folder/degree/flat_min.2.csv"
    run:
        df = pd.read_csv(input.flat_probability,sep="\t")

        bait_degree, params_bait = get_degree(df)
        prey_degree, params_prey = get_degree(df,False)
        full_degree = bait_degree.merge(prey_degree,on="gene_name",how="outer")
        full_degree = fill_na(full_degree,params_bait,params_prey)

        full_degree.to_csv(output.summed_probability,sep="\t",index=False)

        threshold_degree(df,0.1).to_csv(output.threshold_1,sep="\t",index=False)
        threshold_degree(df,0.2).to_csv(output.threshold_2,sep="\t",index=False)


