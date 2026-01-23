import pandas as pd
from itertools import combinations


def subset_ppi(ppi_df, gene_set):
    ppi_df_ss = ppi_df[
        (ppi_df["gene_name_bait"].isin(gene_set)) &
        (ppi_df["gene_name_prey"].isin(gene_set))
        ]
    return ppi_df_ss


def get_genes_in_partitions(partition_df, partition_set):
    genes_partition_set = set()
    for c_partition in partition_set:
        genes_partition_set |= set(partition_df[partition_df["sequence_partition"] == c_partition]["gene_name"])
    return genes_partition_set


def estimate_ppis_kept(ppi_df, partition_df, partition_set):
    all_genes_in_partitions = get_genes_in_partitions(partition_df,partition_set)
    ppis_in_partitions = subset_ppi(ppi_df,all_genes_in_partitions)
    return ppis_in_partitions.shape[0] / ppi_df.shape[0]


def maximise_data_kept(ppi_df, partition_df, remaining_partitions):
    test_slice = int(len(remaining_partitions) * 2 / 3)
    validation_slice = len(remaining_partitions) - test_slice
    test_partitions = list(combinations(remaining_partitions,test_slice))
    combination_partitions = [b for tp in test_partitions for b in remaining_partitions if b not in tp]
    combination_partitions = [
        combination_partitions[i:(i + validation_slice)] for i in range(0,len(combination_partitions),validation_slice)
    ]
    test_partitions, validation_partitions = sorted(
        zip(test_partitions,combination_partitions),
        key=lambda x: estimate_ppis_kept(ppi_df,partition_df,x[0]) + estimate_ppis_kept(ppi_df,partition_df,
            x[1]),reverse=True
    )[0]
    return test_partitions, validation_partitions


rule get_negative_set:
    input:
        input_pod=f"work_folder{pn}/analysis/POD/POD_{{dataset}}.csv",
    output:
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_full_neg_{{neg_limit}}_limit.csv"
    run:
        df_pod = pd.read_csv(input.input_pod,sep="\t")
        df_neg = df_pod[
            (df_pod["n_observed"] == 0) &
            (df_pod["n_tested"] >= int(wildcards.neg_limit))]
        df_neg.to_csv(output.full_neg, sep="\t", index=False)

rule get_positive_set:
    input:
        input_pod=f"work_folder{pn}/analysis/POD/POD_{{dataset}}.csv",
    output:
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_full_limit_poslimit_{{pos_limit}}_pos.csv"
    run:
        df_pod = pd.read_csv(input.input_pod,sep="\t")
        df_pos = df_pod[df_pod["lower_bound_pod"] >= float(wildcards.pos_limit)]
        df_pos.to_csv(output.full_neg, sep="\t", index=False)

rule define_negative_sets:
    input:
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_full_neg_{{neg_limit}}_limit.csv",
        train_partition_genes= f"work_folder{pn}/subsets/train/genes/genes_{{dataset}}_{{pos_limit}}.txt",
        validationq_partition_genes= f"work_folder{pn}/subsets/validation/genes/genes_{{dataset}}_{{pos_limit}}.txt",
        test_partition_genes= f"work_folder{pn}/subsets/test/genes/genes_{{dataset}}_{{pos_limit}}.txt"
    output:
        train_neg = f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        val_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
    run:
        df_neg = pd.read_csv(input.full_neg,sep="\t")

        for gene_file, output_file in zip(
                [input.train_partition_genes, input.validationq_partition_genes, input.test_partition_genes],
                [output.train_neg, output.val_neg, output.test_neg]):
            with open(gene_file, "r") as f:
                genes = {gene.strip() for gene in f}
            subset_ppi(df_neg,genes).to_csv(
                output_file,sep="\t",index=False)


rule define_positive_sets:
    input:
        gene_partition=f"work_folder{pn}/protein_sequences/similarity/{{dataset}}/gene_partition_poslimit_{{pos_limit}}.tsv",
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_full_limit_poslimit_{{pos_limit}}_pos.csv"
    output:
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_limit_{{pos_limit}}_pos.csv",
        train_partition_genes = f"work_folder{pn}/subsets/train/genes/genes_{{dataset}}_{{pos_limit}}.txt",
        val_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_limit_{{pos_limit}}_pos.csv",
        validate_partition_genes= f"work_folder{pn}/subsets/validation/genes/genes_{{dataset}}_{{pos_limit}}.txt",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_limit_{{pos_limit}}_pos.csv",
        test_partition_genes= f"work_folder{pn}/subsets/test/genes/genes_{{dataset}}_{{pos_limit}}.txt",
    run:
        df_pos = pd.read_csv(input.full_pos,sep="\t")
        partitions_df = pd.read_csv(input.gene_partition,sep="\t")
        partitions = partitions_df["sequence_partition"].unique()
        train_combinations = list(combinations(partitions,int(len(partitions) * 0.6)))
        train_partition = sorted(
            train_combinations,key=lambda x: estimate_ppis_kept(df_pos,partitions_df,x),reverse=True
        )[0]

        test_validation_partitions = [p for p in partitions if p not in train_partition]
        validation_partition, test_partition  = maximise_data_kept(df_pos,partitions_df,test_validation_partitions)

        pos_kept_train = estimate_ppis_kept(df_pos,partitions_df,train_partition)
        pos_kept_validate = estimate_ppis_kept(df_pos,partitions_df,validation_partition)
        pos_kept_test = estimate_ppis_kept(df_pos,partitions_df,test_partition)
        msg = f"From the KaFFPa partitions the number of kept {round(100 * (pos_kept_train + pos_kept_test + pos_kept_validate))} % of PPIs: \n"
        msg += f"\tTrain: {round(100 * pos_kept_train)} % \n"
        msg += f"\tValidation: {round(100 * pos_kept_validate)} % \n"
        msg += f"\tTest: {round(100 * pos_kept_test)} % \n"
        print(msg)

        for i, partition in enumerate([train_partition, validation_partition, test_partition]):
            genes = get_genes_in_partitions(partitions_df, partition)
            subset_ppi(df_pos,genes).to_csv(
                output[i * 2],sep="\t",index=False)  # NOTE: If order of output is changed, this breaks
            with open(output[i * 2 + 1], "w") as w:
                _ = [w.write(gene +"\n") for gene in genes]

rule maxflow_splits:
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
        min_max_flow =65
    log: "logs/maxflow/maxflow_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{pos_limit}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    output:
        set_max_flow_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.edgelist",
        set_max_flow_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.edgelist"
    shell:
        """(
        python3 {params.script_location} \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --max_flow_positive {output.set_max_flow_pos} \
            --max_flow_negative {output.set_max_flow_neg} \
            --min_max_flow {params.min_max_flow} \
            --subset {wildcards.settype}
        ) >{log} 2>&1"""

rule ilp:
    params:
        script_location="src/PPIClassification/DataSplit/ILP.py",
        accepted_missmatch = 2
    threads: 20
    log: "logs/ilp/ilp_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_max_flow_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.edgelist",
        set_max_flow_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.edgelist"
    output:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv"
    shell:
        """(
        
        python3 {params.script_location} \
            --positive_data {input.set_max_flow_pos} \
            --negative_data {input.set_max_flow_neg} \
            --balanced_positive {output.balanced_pos} \
            --balanced_negative {output.balanced_neg} \
            --accepted_error {params.accepted_missmatch} \
            --threads {threads} \
            --subset {wildcards.settype}
        ) >{log} 2>&1"""
