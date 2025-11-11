import pandas as pd
from scipy.stats import fisher_exact

rule negatome_comparison:
    """
    Compare if Negatome2.0 is comparable to high confidence non-interactors
    """
    params:
        negatome2="data/PFAM-manual-stringent-negatome2.csv"
    input:
        pod_data=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
        uniprot=f"work_folder{pn}/intact/uniprot_to_gene_name.csv"
    output:
        table=f"work_folder{pn}/analysis/neg2compare/{{data}}.txt"
    run:
        uniprot_2_gene = pd.read_csv(input.uniprot,sep="\t")
        pod_df = pd.read_csv(input.pod_data,sep="\t")
        neg2_df = pd.read_csv(params.negatome2,sep="\t")

        p_cols = ["ProteinA", "ProteinB"]
        neg2_df = neg2_df[p_cols]
        neg2_df_swp = neg2_df.rename(columns={'ProteinA': 'ProteinB', 'ProteinB': 'ProteinA'})[p_cols]

        neg2_directional_df = pd.concat([neg2_df, neg2_df],ignore_index=True).rename(
            columns={
                p: f"uniprot_{c}" for p, c in zip(p_cols,["bait", "prey"])
            })
        neg2_genes = neg2_directional_df.merge(
            uniprot_2_gene,left_on="uniprot_bait",right_on="uniprot_id"
        ).merge(
            uniprot_2_gene,left_on="uniprot_prey",right_on="uniprot_id",suffixes=["_bait", "_prey"]
        )[["gene_name_bait", "gene_name_prey"]]
        neg2_genes["in_neg2"] = True
        pod_neg = pod_df.merge(neg2_genes,on=["gene_name_bait", "gene_name_prey"],how="outer").copy()
        pod_neg.loc[pod_neg["in_neg2"].isna(), "in_neg2"] = False
        pod_neg = pod_neg.loc[~pod_neg["upper_bound_pod"].isna()]
        if wildcards.data == "abundance_mcmc":
            pod_neg["n_tested"] = pod_neg[[c for c in pod_neg.columns if "n_tested_" in c]].sum(axis=1)
            pod_neg["n_observed"] = pod_neg[[c for c in pod_neg.columns if "n_observed_" in c]].sum(axis=1)
        c_table = pod_neg.groupby("in_neg2")[["n_tested", "n_observed"]].sum()
        c_table["n_not_observed"] = c_table["n_tested"] - c_table["n_observed"]
        OR, p_value = fisher_exact(c_table[["n_observed", "n_not_observed"]])
        with open(output.table, "w") as w:
            w.write(str(c_table) + "\n")
            w.write(f"OR: {OR}\n")
            w.write(f"P-value: {p_value}\n")

