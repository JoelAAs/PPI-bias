import pandas as pd


rule define_positive_negative_sets:
    params:
        neg_limit=4,
        pos_limit=0.15 # 2/2 or 3/4 etc
    input:
        gene_partition = f"work_folder{pn}/protein_sequences/similarity/gene_int_id.tsv",
        input_pod = f"work_folder{pn}/analysis/POD/{{dataset}}_ms.csv"
    output:
        train_pos = f"work_folder/{pn}/subsets/train/{{dataset}}}_pos.csv",
        train_neg = f"work_folder/{pn}/subsets/train/{{dataset}}_neg.csv",
        test_pos = f"work_folder/{pn}/subsets/test/{{dataset}}_pos.csv",
        test_neg = f"work_folder/{pn}/subsets/test/{{dataset}}_neg.csv",
        val_pos = f"work_folder/{pn}/subsets/validation/{{dataset}}_pos.csv",
        val_neg = f"work_folder/{pn}/subsets/validation/{{dataset}}_neg.csv"
    run:
        df_tests = pd.read_csv(input.input_pod, sep="\t")
        partitions_df = pd.read_csv(input.gene_partition, sep="\t")

        partitions = partitions_df["partition_id"].unique()

    df_ms = pd.read_csv("work_folder/per_gene/analysis/POD/POD_ms.csv", sep="\t")
    df_neg = df_ms[(df_ms["n_tested"] > 3) & (df_ms["n_observed"] == 0)]
    df_neg = df_neg[["gene_name_bait", "gene_name_prey"]].copy()
    df_neg.columns = ["bait", "prey"]

    df_pos = df_ms[df_ms["lower_bound_pod"] > 0.2]
    df_pos = df_pos[["gene_name_bait", "gene_name_prey"]].copy()
    df_pos.columns = ["bait", "prey"]

    subset_negative_set(df_neg, df_pos)
