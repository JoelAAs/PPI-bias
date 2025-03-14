import gzip
import pandas as pd
import numpy as np
from formating import *


###
def get_index_dict(filename, dtypes = (str, str)):
    index_dict = dict()
    with open(filename, "r") as f:
        for l in f:
            key, value = l.strip().split("\t")
            index_dict[dtypes[0](key)] = dtypes[1](value)

    return index_dict

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

rule get_protein_index:
    input:
        formated = "work_folder/intact/method_subset/bait_prey_publications.csv"
    output:
        protein_to_index_file = "work_folder/intact/index/protein_to_index.csv",
        index_to_protein_file = "work_folder/intact/index/index_to_protein.csv"
    run:
        bait_prey_df = pd.read_csv(input.formated, sep = "\t")
        all_proteins = sorted(list(set(
            bait_prey_df["bait"].tolist() +
            bait_prey_df["prey"].tolist()
        )))
        with open(output.index_to_protein_file, "w") as wip:
            with open(output.protein_to_index_file, "w") as wpi:
                for i, protein in enumerate(all_proteins):
                    wip.write(f"{i}\t{protein}\n")
                    wpi.write(f"{protein}\t{i}\n")

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
    input:
        pubmed_id = "work_folder/pid_cell_line/{cell_line}.csv",
        method_subset = "work_folder/intact/method_subset/{method}.csv"
    output:
        bait_prey_table = "work_folder/intact/method_subset/cell_line_subset/interactions_{cell_line}_{method}.csv",
    run:
        with open(input.pubmed_id, "r") as f:
            pids = [l.strip() for l in f]
            pids = pids[1:]

        method_df = pd.read_csv(input.method_subset, sep = "\t")
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
        pod_file = config["cell_lines"][wildcards.cell_line]["pod"]

        max_interaction_dict, observed_interaction_dict  = get_interaction_dict(
            bait_prey_df,
            method=method,
            prey_file=pod_file)

        dict_to_pairs_file(
            max_interaction_dict,
            observed_interaction_dict,
            output.protein_pairs,
            method = method,
            prey_pod_file = pod_file)
