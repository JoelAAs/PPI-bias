from collections import defaultdict

def nested_dict():
    return defaultdict(nested_dict)

def set_or_update_nested(c_dict, keys, value):
    for key in keys[:-1]:
        c_dict = c_dict[key]
    if keys[-1] in c_dict:
        c_dict[keys[-1]] += value
    else:
        c_dict[keys[-1]] = value

def write_dict_to_file(current_dict, file_w, line):
    if "n_tested" in current_dict:
        file_w.write(
            line + f"{current_dict["n_tested"]}\t{current_dict["n_observed"]}\n")
    else:
        for key, next_dict in current_dict.items():
            write_dict_to_file(next_dict, file_w, line + f"{key}\t")

def aggregate_inferred_experiments(pids, output_file, cl=False):
    """
    Aggregated counts of tests per experiment
    :param pids: (list) of files to be aggregated
    :param output_file: (str) of output file location
    :param cl: (boolean) if cell_line should be considered or not
    :return: -
    """
    ppi_dict = nested_dict()
    for study in pids:
        with open(study, "r") as f:
            header = True
            for line in f:
                if header:
                    header = False
                else:
                    values = line.strip().split("\t")
                    bait, prey, n_tested, n_observed = values[:4]
                    if cl:
                        cl_id = values[-1]
                        keys = [bait, prey, cl_id]
                    else:
                        keys = [bait, prey]

                    set_or_update_nested(ppi_dict, keys +["n_tested"], int(n_tested))
                    set_or_update_nested(ppi_dict, keys +["n_observed"], int(n_observed))

    with open(output_file, "w") as w:
        if cl:
            header_line = "gene_name_bait\tgene_name_prey\tcl_id\tn_tested\tn_observed\n"
        else:
            header_line = "gene_name_bait\tgene_name_prey\tn_tested\tn_observed\n"
        w.write(header_line)
        write_dict_to_file(ppi_dict, w, "")