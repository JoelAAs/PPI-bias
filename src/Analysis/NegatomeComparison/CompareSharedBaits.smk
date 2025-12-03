from collections import Counter
import pandas as pd
import networkx as nx
from scipy.stats import entropy

rule get_bait_occurrence:
    input:
        pod=""
    output:
        edge_list_shared=f"work_folder{pn}/analysis/negatome/bait_shared_studies_{{data}}.csv"
    run:
        df_pod = pd.read_csv(input.pod,sep="\t")
        bait_usage = df_pod[~df_pod[["gene_name_bait", "pubmed_id"]].duplicated(keep="first")].copy()
        bait_usage["pubmed_id"] = bait_usage["pubmed_id"].apply(lambda x: set(x.split(";")))

        bait_usage = bait_usage.groupby("gene_name_bait")["pubmed_id"].apply(lambda x: set.union(*x))
        bait_pub_list = list(bait_usage.items())

        with open(output.edge_list_shared,"w") as w:
            w.write("bait_a\tbait_b\tn_shared\n")
            for i, (bait_a, studies_a) in enumerate(bait_pub_list):
                j = i + 1
                if j != len(bait_pub_list):
                    for bait_b, studies_b in bait_pub_list[j:]:
                        n_shared_studies = len(studies_a & studies_b)
                        w.write(f"{bait_a}\t{bait_b}\t{n_shared_studies}\n")


rule get_bait_bait_degree:
    params:
        limits=[2, 3, 4, 5]
    input:
        edge_list_shared=f"work_folder{pn}/analysis/negatome/bait_shared_studies_{{data}}.csv"
    output:
        bait_bait_degree=f"work_folder{pn}/analysis/negatome/bait_bait_degree_{{data}}.csv"
    run:
        df_bait_bait = pd.read_csv(input.edge_list_shared,sep="\t")

        limits_degree = []
        for limit in params.limits:
            ss_edges = df_bait_bait[df_bait_bait["n_shared"] > limit]
            G = nx.from_pandas_edgelist(ss_edges,source="bait_a",target="bait_b")
            limits_degree.append(
                pd.DataFrame(dict(G.degree()).items(),columns=["gene", f"bait_degree_{limit}"]).set_index("gene")
            )
        joined_df = limits_degree[0].join(limits_degree[1:],how='outer').reset_index().fillna(0)
        joined_df.to_csv(output.bait_bait_degree,sep="\t",index=False)


rule negative_data_vs_bait_degree:
    params:
        limits=[2, 3, 4, 5]
    input:
        pod=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
        bait_bait_degree=f"work_folder{pn}/analysis/negatome/bait_bait_degree_{{data}}.csv"
    output:
        neg_bait_degree=f"work_folder{pn}/analysis/negatome/neg_bait_bait_degree_{{data}}.csv"
    run:
        df_pod = pd.read_csv(input.pod,sep="\t")
        df_pod_neg = df_pod[df_pod["n_observed"] == 0]
        bait_bait_degree = pd.read_csv(input.bait_bait_degree,sep="\t")
        negative_degree = []
        cols = ["gene"]
        for limit in params.limits:
            ss_df_pod_neg = df_pod_neg[df_pod_neg["n_tested"] > limit]
            ss_count = ss_df_pod_neg.groupby("gene_name_bait").size()
            negative_degree.append(ss_count)
            cols.append(f"negative_degree_{limit}")

        negative_degree_df = pd.concat(negative_degree,axis=1).fillna(0)
        negative_degree_df = negative_degree_df.reset_index()
        negative_degree_df.columns = cols

        bait_negative_degree = negative_degree_df.merge(bait_bait_degree,on="gene",how="outer").fillna(0)
        bait_negative_degree.to_csv(output.neg_bait_degree,sep="\t",index=False)


rule non_interaction_prey_entropy_entropy:
    input:
        pod=f"work_folder{pn}/analysis/POD/POD_{{data}}.csv",
        bait_prey_degree=f"work_folder{pn}/formated/bait_prey_publications.csv"
    output:
        entropy_annotated=f"work_folder{pn}/analysis/negatome/test_entropy_{{data}}_limit_{{min_tests}}.csv"
    run:
        df_pod = pd.read_csv(input.pod,sep="\t")
        df_pod_neg = df_pod[df_pod["n_observed"] == 0]
        df_pod_neg = df_pod_neg[df_pod_neg["n_tested"] >= int(wildcards.min_tests)]

        bait_prey_df = pd.read_csv(input.bait_prey_degree,sep="\t")
        bait_prey_df["pub_method"] = bait_prey_df.apply(
            lambda x: "_".join([str(x["pubmed_id"]), x["detection_method"]]),axis=1)

        write_header = True
        for pids_all in df_pod_neg["pubmed_id"].unique():
            protein_pair_pub_ss = df_pod_neg[df_pod_neg["pubmed_id"] == pids_all]
            pids = pids_all.split(";")
            bait_prey_studies = bait_prey_df[bait_prey_df["pub_method"].isin(pids)]

            for prey in protein_pair_pub_ss["gene_name_prey"].unique():
                prey_bait_ss = bait_prey_studies[
                    bait_prey_studies["gene_name_prey"] == prey][["gene_name_bait", "pub_method"]]
                study_normilised_count = prey_bait_ss.groupby(
                    "pub_method")["gene_name_bait"].value_counts(normalize=True).groupby("gene_name_bait").sum()

                prey_pub_entropy = entropy(study_normilised_count)
                protein_pair_prey_ss = protein_pair_pub_ss[protein_pair_pub_ss["gene_name_prey"] == prey].copy()

                protein_pair_prey_ss["pair_entropy"] = prey_pub_entropy
                if write_header:
                    protein_pair_prey_ss.to_csv(output.entropy_annotated,sep="\t",index=False)
                    write_header = False
                else:
                    protein_pair_prey_ss.to_csv(output.entropy_annotated,sep="\t",mode="a",index=False,header=False)
