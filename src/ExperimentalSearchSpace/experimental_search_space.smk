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
            create_or_update(tested_bait_prey_dict[bait],prey,1) for prey in preys
        ]

    for _, (bait, prey) in ppi_ss[[f"{id_pattern}_bait", f"{id_pattern}_prey"]].iterrows():
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


checkpoint infer_experimental_search_space:
    params:
        id_pattern = config["id_pattern"],
        cell_line_present = config["cell_line_present"]
    input:
        bait_prey_file = config["ppi_df"]
    output:
        directory("work_folder/inferred_search_space/experimental")
    run:
        os.makedirs(output[0], exist_ok=True)
        bait_prey_df = pd.read_csv(input.bait_prey_file, sep="\t")
        bait_prey_df = bait_prey_df[
            bait_prey_df[f"{params.id_pattern}_bait"] != bait_prey_df[f"{params.id_pattern}_prey"]
        ]
        for pid in bait_prey_df["pubmed_id"].unique():
            pid_ss = bait_prey_df[bait_prey_df["pubmed_id"] == pid]
            for detection_method in pid_ss["detection_method"].unique():
                tested_bait_prey_dict      = dict()
                observation_bait_prey_dict = dict()
                method_pid_ss = pid_ss[
                    pid_ss["detection_method"] == detection_method
                ]
                id_cols = [f"{params.id_pattern}_baits", f"{params.id_pattern}_prey"]
                method_pid_ss = method_pid_ss[~method_pid_ss[id_cols].duplicated()] # Isoforms iof gene name gives more observed than tested
                if params.cell_line_present:
                    for cl_id in method_pid_ss["cl_id"].unique():
                        output_file = f"work_folder/inferred_search_space/experimental/{pid}_{detection_method}_{cl_id}.csv"
                        cl_method_pid_ss = method_pid_ss[method_pid_ss["cl_id"] == cl_id]
                        tested_dict, obs_dict = get_tested_observed_dicts(cl_method_pid_ss, params.id_pattern)
                        write_observed(
                            tested_dict=tested_dict, obs_dict=obs_dict,
                            output_file=output_file, id_pattern=params.id_pattern,
                            detection_method=detection_method, pid=pid, cl=cl_id)


                else:
                    output_file = f"work_folder/inferred_search_space/experimental/{pid}_{detection_method}.csv"
                    tested_dict, obs_dict = get_tested_observed_dicts(method_pid_ss, params.id_pattern)
                    write_observed(
                        tested_dict=tested_dict,obs_dict=obs_dict,
                        output_file=output_file,id_pattern=params.id_pattern,
                        detection_method=detection_method,pid=pid)




