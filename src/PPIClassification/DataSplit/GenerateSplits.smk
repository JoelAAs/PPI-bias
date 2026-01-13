import pandas as pd


rule define_positive_negative_sets:
    params:
        neg_limit=4,
        pos_limit=0.15 # 2/2 or 3/4 etc
    input:
        gene_partition = f"work_folder{pn}/protein_sequences/similarity/gene_partition.tsv",
        input_pod = f"work_folder{pn}/analysis/POD/POD_{{dataset}}.csv"
    output:
        train_pos = f"work_folder{pn}/subsets/train/{{dataset}}_pos.csv",
        train_neg = f"work_folder{pn}/subsets/train/{{dataset}}_neg.csv",
        test_pos = f"work_folder{pn}/subsets/test/{{dataset}}_pos.csv",
        test_neg = f"work_folder{pn}/subsets/test/{{dataset}}_neg.csv",
        val_pos = f"work_folder{pn}/subsets/validation/{{dataset}}_pos.csv",
        val_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_neg.csv"
    run:
        df_tests = pd.read_csv(input.input_pod, sep="\t")
        partitions_df = pd.read_csv(input.gene_partition, sep="\t")
        df_pos = df_tests[df_tests["lower_bound_pod"] >= params.pos_limit]
        df_neg = df_tests[(df_tests["n_observed"] == 0) & (df_tests["n_tested"] >= params.neg_limit)]
        partitions = partitions_df["sequence_partition"].unique()
        partition_pos = dict()
        partition_neg = dict()
        for partition_id in partitions:
            genes_partition = set(partitions_df[partitions_df["sequence_partition"] == partition_id]["gene_name"])
            partition_pos[partition_id] = df_pos[
                (df_pos["gene_name_bait"].isin(genes_partition)) & (df_pos["gene_name_prey"].isin(genes_partition))
            ]
            partition_neg[partition_id] = df_neg[
                (df_neg["gene_name_bait"].isin(genes_partition)) &( df_neg["gene_name_prey"].isin(genes_partition))
                ]
        # Largest to train, test, validate
        order_partitions = sorted(partitions, key=lambda x: partition_pos[x].shape[0], reverse=True) # Will break if k!=3 kahip
        for i, partition_id in enumerate(order_partitions):
            partition_pos[partition_id].to_csv(output[i*2], sep="\t", index=False)
            partition_neg[partition_id].to_csv(output[i*2+1],sep="\t",index=False)

rule balance_splits:
    input:
        set_pos = f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_pos.csv",
        set_neg= f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_neg.csv"
    output:
        set_balanced_pos = f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_pos.csv",
        set_balanced_neg = f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_neg.csv"
    shell:
        """
        
        """