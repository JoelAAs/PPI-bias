import pandas as pd
import networkx as nx

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
        neg_bait_degree = f"work_folder{pn}/analysis/negatome/neg_bait_bait_degree_{{data}}.csv"
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

        negative_degree_df = pd.concat(negative_degree, axis = 1).fillna(0)
        negative_degree_df = negative_degree_df.reset_index()
        negative_degree_df.columns = cols

        bait_negative_degree = negative_degree_df.merge(bait_bait_degree, on="gene", how="outer").fillna(0)
        bait_negative_degree.to_csv(output.neg_bait_degree, sep="\t", index=False)