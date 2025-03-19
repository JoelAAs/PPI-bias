import gzip
from fileinput import input

import pandas as pd
import numpy as np
from formating import *



##### Rules
rule format_intact:
    input:
        intact = "data/intact/human.txt"
    output:
        formated = "work_folder/intact/method_subset/bait_prey_publications.csv"
    run:
        intact_df = read_intact(input.intact)
        bait_prey_df = filter_and_index(intact_df)
        bait_prey_df.to_csv(
            output.formated,
            sep="\t",
            index = None
        )

rule get_gene_name_uniprot:
    input:
        intact = "data/intact/human.txt"
    output:
        uniprot = "work_folder/intact/uniprot_to_gene_name.csv"
    run:
        get_gene_names(input.intact, output.uniprot)

rule subset_method:
    input:
        formated = "work_folder/intact/method_subset/bait_prey_publications.csv"
    output:
        expand(
            "work_folder/intact/method_subset/{method}.csv", method=config["methods"]
        )
    run:
        intact_df = pd.read_csv(input.formated, sep="\t")
        for method in config["methods"]:
            method_subset = intact_df[intact_df["detection_method"] == method]
            method_subset.to_csv(
                f"work_folder/intact/method_subset/{method}.csv",
                sep="\t",
                index=False
            )

rule subset_cell_line:
    params:
        max_cl_mentions = 2
    input:
        pubmed_id = "work_folder/pid_cell_line/{cell_line}.csv",
        method_subset = "work_folder/intact/method_subset/{method}.csv"
    output:
        bait_prey_table = "work_folder/intact/method_subset/cell_line_subset/interactions_{cell_line}_{method}.csv",
    run:
            pid_df = pd.read_csv(input.pubmed_id,
                sep = "\t",
                dtype={"pubmed_id":str, "cl_count":int}
            )
            pid_df = pid_df[pid_df["cl_count"] <= params.max_cl_mentions]
            pids = pid_df["pubmed_id"].unique()
            method_df = pd.read_csv(
                input.method_subset,
                sep = "\t",
                dtype="str"
            )
            method_df = method_df[method_df["pubmed_id"].isin(pids)]
            method_df.to_csv(output.bait_prey_table, sep = "\t", index = False)


rule get_ppi_counts:
    input:
        bait_prey_table = "work_folder/intact/method_subset/cell_line_subset/interactions_{cell_line}_{method}.csv",
        prey_pod = "work_folder/cell_type_pod/{cell_line}_pod.csv"
    output:
        protein_pairs = "work_folder/intact/pair_count/ppi_pair_counts_{cell_line}_{method}.csv"

    run:
        bait_prey_df = pd.read_csv(input.bait_prey_table, sep="\t")
        method =  config["methods"][wildcards.method]

        max_interaction_dict, observed_interaction_dict  = get_interaction_dict(
            bait_prey_df=bait_prey_df,
            method=method,
            prey_file=input.prey_pod)

        dict_to_pairs_file(
            max_interaction_dict=max_interaction_dict,
            observed_interaction_dict=observed_interaction_dict,
            output_filename=output.protein_pairs,
            method = method,
            prey_pod_file = input.prey_pod)
