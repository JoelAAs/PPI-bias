import pandas as pd
from itertools import combinations

def get_gene_partition(wc):
    if wc["partition_name"] == "sequencesimilarity":
        return f"work_folder{pn}/subsets/partitions/sequencesimilarity_gene_name.txt"
    elif wc["partition_name"] == "maxpos":
        return f"work_folder{pn}/subsets/partitions/{wc['dataset']}_{wc["network_type"]}_limit_{wc['pos_limit']}_gene_name.txt"
    else:
        raise ValueError(f"unknown partition name = {wc['partition_name']}")


def subset_ppi(ppi_df, gene_set):
    ppi_df_ss = ppi_df[
        (ppi_df["gene_name_bait"].isin(gene_set)) &
        (ppi_df["gene_name_prey"].isin(gene_set))
        ]
    return ppi_df_ss


def get_genes_in_partitions(partition_df, partition_set):
    genes_partition_set = set()
    for c_partition in partition_set:
        genes_partition_set |= set(partition_df[partition_df["partition"] == c_partition]["gene_name"])
    return genes_partition_set


def estimate_ppis_kept(ppi_df, partition_df, partition_set):
    all_genes_in_partitions = get_genes_in_partitions(partition_df,partition_set)
    ppis_in_partitions = subset_ppi(ppi_df,all_genes_in_partitions)
    return ppis_in_partitions.shape[0] / ppi_df.shape[0]


def maximise_data_kept(ppi_df, partition_df, remaining_partitions):

    validation_slice = int(len(remaining_partitions) * 1/2)
    test_slice = len(remaining_partitions) - validation_slice
    validation_partitions = list(combinations(remaining_partitions,validation_slice))
    combination_partitions = [rp for vp in validation_partitions for rp in remaining_partitions if rp not in vp]
    combination_partitions = [
        combination_partitions[i:(i + validation_slice)] for i in range(0,len(combination_partitions),test_slice)
    ]
    test_partitions, validation_partitions = sorted(
        zip(validation_partitions, combination_partitions),
        key=lambda x: estimate_ppis_kept(
            ppi_df,partition_df,x[0]) + estimate_ppis_kept(ppi_df,partition_df,x[1]),reverse=True
    )[0]
    return test_partitions, validation_partitions


rule get_negative_set:
    input:
        input_pod=f"work_folder{pn}/analysis/POD/{{network_type}}/POD_{{dataset}}.pq",
    output:
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_full_{{neg_limit}}_neg.pq"
    resources:
        mem_gb=50
    run:
        df_pod = pd.read_parquet(input.input_pod)
        df_neg = df_pod[
            (df_pod["n_observed"] == 0) &
            (df_pod["n_tested"] >= int(wildcards.neg_limit))]
        df_neg.to_parquet(output.full_neg, index=False)

rule get_positive_set:
    input:
        input_pod=f"work_folder{pn}/analysis/POD/{{network_type}}/POD_{{dataset}}.pq",
    output:
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_full_{{pos_limit}}_pos.pq"
    resources:
        mem_gb=50
    run:
        df_pod = pd.read_parquet(input.input_pod)
        df_pos = df_pod[df_pod["lower_bound_pod"] >= float(wildcards.pos_limit)]
        df_pos.to_parquet(output.full_pos, index=False)


rule define_partitions:
    input:
        gene_partition = lambda wc: get_gene_partition(wc),
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_full_{{pos_limit}}_pos.pq"
    output:
        train_partition_genes = f"work_folder{pn}/subsets/train/genes/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        validate_partition_genes= f"work_folder{pn}/subsets/validation/genes/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        test_partition_genes= f"work_folder{pn}/subsets/test/genes/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
    run:
        df_pos = pd.read_parquet(input.full_pos)
        partitions_df = pd.read_csv(input.gene_partition, sep="\t")
        partitions = partitions_df["partition"].unique()
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
            with open(output[i], "w") as w:
                _ = [w.write(gene +"\n") for gene in genes]


rule define_positive_split:
    input:
        full_pos = f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_full_{{pos_limit}}_pos.pq",
        train_partition_genes = f"work_folder{pn}/subsets/train/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        validate_partition_genes= f"work_folder{pn}/subsets/validation/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        test_partition_genes= f"work_folder{pn}/subsets/test/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt"
    output:
        train_pos = f"work_folder{pn}/subsets/train/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_{{partition_name}}_pos.pq",
        val_pos = f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_{{partition_name}}_pos.pq",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_{{partition_name}}_pos.pq"
    run:
        df_pos = pd.read_parquet(input.full_pos)
        for partition_file, output_file in zip(
            [input.train_partition_genes, input.validate_partition_genes, input.test_partition_genes],
            [output.train_pos, output.val_pos, output.test_pos]):
            with open(partition_file, "r") as f:
                genes = {gene.strip() for gene in f}
            subset_ppi(df_pos,genes).to_parquet(
                output_file,index=False)



rule define_negative_sets:
    input:
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_{{network_type}}_full_{{neg_limit}}_neg.pq",
        train_partition_genes= f"work_folder{pn}/subsets/train/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        validation_partition_genes= f"work_folder{pn}/subsets/validation/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt",
        test_partition_genes= f"work_folder{pn}/subsets/test/genes/cdhit/genes_{{dataset}}_{{network_type}}_{{pos_limit}}_{{partition_name}}.txt"
    output:
        train_neg = f"work_folder{pn}/subsets/train/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq",
        val_neg = f"work_folder{pn}/subsets/validation/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq"
    resources:
        mem_gb=80
    run:
        df_neg = pd.read_parquet(input.full_neg)
        print("Rows:", len(df_neg))
        print("Memory GB:", df_neg.memory_usage(deep=True).sum() / 1e9, flush=True)
        genes_partitions = [
            {gene.strip() for gene in open(f)} for f in [
                input.train_partition_genes, input.validation_partition_genes, input.test_partition_genes
                ]
            ]
        df_neg["gene_name_bait"] = df_neg["gene_name_bait"].astype("category")
        df_neg["gene_name_prey"] = df_neg["gene_name_prey"].astype("category")
        baits = df_neg["gene_name_bait"]
        prey = df_neg["gene_name_prey"]


        for set_partition, output_file in zip(
                gene_partitions
                [output.train_neg, output.val_neg, output.test_neg]):

            partition_mask = (
                (baits.isin(set_partition)) &
                (prey.isin(set_partition))
            )
            df_neg[partition_mask].to_parquet(
                output_file,index=False)