import numpy as np
import pandas as pd
import os


def create_or_update(c_dict, key, value):
    if key in c_dict:
        if value != dict():
            c_dict[key] += value
    else:
        c_dict[key] = value

def get_tested_observed_dicts(ppi_ss, id_pattern, include_self_interactions=False):
    tested_bait_prey_dict = dict()
    observation_bait_prey_dict = dict()

    preys = ppi_ss[f"{id_pattern}_prey"].unique()
    for bait in ppi_ss[f"{id_pattern}_bait"].unique():
        create_or_update(tested_bait_prey_dict,bait,dict())

        _ = [
            create_or_update(tested_bait_prey_dict[bait],prey,1)
            for prey in preys if prey != bait or include_self_interactions
        ]

    for _, (bait, prey) in ppi_ss[[f"{id_pattern}_bait", f"{id_pattern}_prey"]].iterrows():
        if bait != prey or include_self_interactions:
            create_or_update(observation_bait_prey_dict,bait,dict())
            create_or_update(observation_bait_prey_dict[bait],prey,1)

    return tested_bait_prey_dict, observation_bait_prey_dict

def write_observed(tested_dict, obs_dict, output_file, id_pattern, detection_method, pid, cl=np.nan):
    header = "\t".join([
        f"{id_pattern}_bait",
        f"{id_pattern}_prey",
        "n_tested",
        "n_observed",
        "detection_method",
        "pubmed_id",
        "CVCL"
    ]) + "\n"
    with open(output_file,"w") as w:
        w.write(header)
        for bait in tested_dict:
            for prey in tested_dict[bait]:
                try:
                    observed = obs_dict[bait][prey]
                except KeyError:
                    observed = 0

                n_tests = tested_dict[bait][prey]

                lineout = (
                    f"{bait}\t{prey}\t{n_tests}\t{observed}\t{detection_method}\t{pid}\t{cl}\n")
                w.write(lineout)

def get_input_ppi_file(cell_line_wc):
    if cell_line_wc == "_cell_line":
        return "work_folder/formated/bait_prey_CVCL.csv"

    else:
        return "work_folder/formated/bait_prey_publications.csv"


checkpoint infer_experimental_search_space:
    """
    Inferres negative data for baits of other preys seen in multibait experiments.
    If {cell_line} is "cell_line" it will infere using cell line specific data specified in config "cell_line_ppis"
    """
    params:
        id_pattern = config["id_pattern"],
        include_y2h_self_interactions = config["include_y2h_self_interactions"],
        y2h_methods = config["y2h"],
        drop_isoforms = config["drop_isoforms"]
    input:
        bait_prey_file = lambda wc: storage.fs(get_input_ppi_file(wc.cell_line))
    output:
        directory("work_folder/inferred_search_space/experimental{cell_line}")
    log:
        "logs/inferred_search_space/experimental{cell_line}.log"
    run:
        os.makedirs(output[0], exist_ok=True)
        bait_prey_df = pd.read_csv(input.bait_prey_file, sep="\t")

        if params.id_pattern == "gene_name":
            id_cols = [
                f"{params.id_pattern}_bait", f"{params.id_pattern}_prey",
                "pubmed_id", "detection_method"
            ]
            if wildcards.cell_line == "_cell_line":
                id_cols.append("CVCL")
            if params.drop_isoforms:
                if params.id_pattern == "uniprot_id": # Drop all isoform info
                    bait_prey_df[f"{params.id_pattern}_bait"] = bait_prey_df[f"{params.id_pattern}_bait"].str.split("-").str[0]
                    bait_prey_df[f"{params.id_pattern}_prey"] = bait_prey_df[f"{params.id_pattern}_prey"].str.split("-").str[0]

            bait_prey_df = bait_prey_df[
                ~bait_prey_df[id_cols].duplicated(keep="first")]  # Isoforms of gene name gives more observed than tested

        for pid in bait_prey_df["pubmed_id"].unique():
            pid_ss = bait_prey_df[bait_prey_df["pubmed_id"] == pid]

            for detection_method in pid_ss["detection_method"].unique():
                method_pid_ss = pid_ss[
                    pid_ss["detection_method"] == detection_method
                ]
                include_current_self_interactions = (
                    detection_method in params.y2h_methods and params.include_y2h_self_interactions
                )

                if wildcards.cell_line == "_cell_line":
                    for cl_id in method_pid_ss["CVCL"].unique():
                        output_file = f"{output[0]}/{pid}_{detection_method}_{cl_id}.csv"
                        cl_method_pid_ss = method_pid_ss[method_pid_ss["CVCL"] == cl_id]
                        tested_dict, obs_dict = get_tested_observed_dicts(
                            cl_method_pid_ss, params.id_pattern, include_current_self_interactions)
                        write_observed(
                            tested_dict=tested_dict, obs_dict=obs_dict,
                            output_file=output_file, id_pattern=params.id_pattern,
                            detection_method=detection_method, pid=pid, cl=cl_id)


                else:
                    output_file = f"{output[0]}/{pid}_{detection_method}.csv"
                    tested_dict, obs_dict = get_tested_observed_dicts(
                        method_pid_ss, params.id_pattern, include_current_self_interactions)
                    write_observed(
                        tested_dict=tested_dict,obs_dict=obs_dict,
                        output_file=output_file,id_pattern=params.id_pattern,
                        detection_method=detection_method,pid=pid)




