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


rule define_positive_negative_sets:
    params:
        neg_limit=2,
        pos_limit=0.15  # 2/2 or 3/4 etc
    input:
        gene_partition=f"work_folder{pn}/protein_sequences/similarity/gene_partition.tsv",
        input_pod=f"work_folder{pn}/analysis/POD/POD_{{dataset}}.csv"
    output:
        train_pos=f"work_folder{pn}/subsets/train/{{dataset}}_pos.csv",
        train_neg=f"work_folder{pn}/subsets/train/{{dataset}}_neg.csv",
        test_pos=f"work_folder{pn}/subsets/test/{{dataset}}_pos.csv",
        test_neg=f"work_folder{pn}/subsets/test/{{dataset}}_neg.csv",
        val_pos=f"work_folder{pn}/subsets/validation/{{dataset}}_pos.csv",
        val_neg=f"work_folder{pn}/subsets/validation/{{dataset}}_neg.csv",
        full_pos=f"work_folder{pn}/subsets/{{dataset}}_full_pos.csv",
        full_neg=f"work_folder{pn}/subsets/{{dataset}}_full_neg.csv"
    run:
        df_tests = pd.read_csv(input.input_pod,sep="\t")
        partitions_df = pd.read_csv(input.gene_partition,sep="\t")
        df_pos = df_tests[df_tests["lower_bound_pod"] >= params.pos_limit]
        df_pos.to_csv(output.full_pos,sep="\t",index=False)
        df_neg = df_tests[(df_tests["n_observed"] == 0) & (df_tests["n_tested"] >= params.neg_limit)]
        df_neg.to_csv(output.full_neg,sep="\t",index=False)
        partitions = partitions_df["sequence_partition"].unique()
        train_combinations = list(combinations(partitions,int(len(partitions) * 0.6)))
        train_partition = sorted(
            train_combinations,key=lambda x: estimate_ppis_kept(df_pos,partitions_df,x),reverse=True
        )[0]

        test_validation_partitions = [p for p in partitions if p not in train_partition]
        test_partition, validation_partition = maximise_data_kept(df_pos,partitions_df,test_validation_partitions)

        pos_kept_train = estimate_ppis_kept(df_pos,partitions_df,train_partition)
        pos_kept_test = estimate_ppis_kept(df_pos,partitions_df,test_partition)
        pos_kept_validate = estimate_ppis_kept(df_pos,partitions_df,validation_partition)
        msg = f"From the KaFFPa partitions the number of kept {round(100 * (pos_kept_train + pos_kept_test + pos_kept_validate))} % of PPIs: \n"
        msg += f"\tTrain: {round(100 * pos_kept_train)} % \n"
        msg += f"\tTest: {round(100 * pos_kept_test)} % \n"
        msg += f"\tValidation: {round(100 * pos_kept_validate)} % \n"
        print(msg)

        for i, partition in enumerate([train_partition, test_partition, validation_partition]):
            genes = get_genes_in_partitions(partitions_df,partition)
            subset_ppi(df_pos,genes).to_csv(
                output[i * 2],sep="\t",index=False)  # NOTE: If order of output is changed, this breaks
            subset_ppi(df_neg,genes).to_csv(output[i * 2 + 1],sep="\t",index=False)

rule maxflow_splits:
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
        min_max_flow =65
    log: "logs/maxflow/maxflow_{settype}_{dataset}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_neg.csv"
    output:
        set_max_flow_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_pos.edgelist",
        set_max_flow_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_neg.edgelist"
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
    log: "logs/ilp/ilp_{settype}_{dataset}.log"
    input:
        set_max_flow_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_pos.edgelist",
        set_max_flow_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_neg.edgelist"
    output:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_neg.csv"
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
