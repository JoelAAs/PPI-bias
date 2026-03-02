import numpy as np
import pandas as pd
import os


def create_or_update(c_dict, key, value):
    if key in c_dict:
        if value != dict():
            c_dict[key] += value
    else:
        c_dict[key] = value

def get_tested_observed_dicts(ppi_ss, id_pattern):
    tested_bait_prey_dict = dict()
    observation_bait_prey_dict = dict()

    preys = ppi_ss[f"{id_pattern}_prey"].unique()
    for bait in ppi_ss[f"{id_pattern}_bait"].unique():
        create_or_update(tested_bait_prey_dict,bait,dict())

        _ = [
            create_or_update(tested_bait_prey_dict[bait],prey,1) for prey in preys if prey != bait
        ]

    for _, (bait, prey) in ppi_ss[[f"{id_pattern}_bait", f"{id_pattern}_prey"]].iterrows():
        if bait != prey:
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
        "cl_id"
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
        return f"work_folder{pn}/formated/bait_prey_CVCL.csv"

    else:
        return f"work_folder{pn}/formated/bait_prey_publications.csv"


checkpoint infer_experimental_search_space:
    """
    Inferres negative data for baits of other preys seen in multibait experiments.
    If {cell_line} is "cell_line" it will infere using cell line specific data specified in config "cell_line_ppis" 
    """
    params:
        id_pattern = config["id_pattern"]
    input:
        bait_prey_file = lambda wc: get_input_ppi_file(wc.cell_line)
    output:
        directory(f"work_folder{pn}/inferred_search_space/experimental{{cell_line}}")
    run:
        os.makedirs(output[0], exist_ok=True)
        bait_prey_df = pd.read_csv(input.bait_prey_file, sep="\t")
        bait_prey_df = bait_prey_df[
            bait_prey_df[f"{params.id_pattern}_bait"] != bait_prey_df[f"{params.id_pattern}_prey"]
        ] 
        if params.id_pattern == "gene_name":
            id_cols = [
                f"{params.id_pattern}_bait", f"{params.id_pattern}_prey",
                "pubmed_id", "detection_method"
            ]
            if wildcards.cell_line == "_cell_line":
                id_cols.append("cl_id")
            bait_prey_df = bait_prey_df[
                ~bait_prey_df[id_cols].duplicated(keep="first")]  # Isoforms iof gene name gives more observed than tested

        for pid in bait_prey_df["pubmed_id"].unique():
            pid_ss = bait_prey_df[bait_prey_df["pubmed_id"] == pid]
            for detection_method in pid_ss["detection_method"].unique():
                tested_bait_prey_dict      = dict()
                observation_bait_prey_dict = dict()
                method_pid_ss = pid_ss[
                    pid_ss["detection_method"] == detection_method
                ]

                if wildcards.cell_line == "_cell_line":
                    for cl_id in method_pid_ss["cl_id"].unique():
                        output_file = f"{output[0]}/{pid}_{detection_method}_{cl_id}.csv"
                        cl_method_pid_ss = method_pid_ss[method_pid_ss["cl_id"] == cl_id]
                        tested_dict, obs_dict = get_tested_observed_dicts(cl_method_pid_ss, params.id_pattern)
                        write_observed(
                            tested_dict=tested_dict, obs_dict=obs_dict,
                            output_file=output_file, id_pattern=params.id_pattern,
                            detection_method=detection_method, pid=pid, cl=cl_id)


                else:
                    output_file = f"{output[0]}/{pid}_{detection_method}.csv"
                    tested_dict, obs_dict = get_tested_observed_dicts(method_pid_ss, params.id_pattern)
                    write_observed(
                        tested_dict=tested_dict,obs_dict=obs_dict,
                        output_file=output_file,id_pattern=params.id_pattern,
                        detection_method=detection_method,pid=pid)




