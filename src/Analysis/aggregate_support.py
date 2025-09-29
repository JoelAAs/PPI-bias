from collections import defaultdict


def update_idx(c_dict, key):
    if key not in c_dict:
        c_dict[key] = len(c_dict) + 1
    return c_dict[key]


def aggregate_inferred_experiments(pids, output_file, cl=False):
    """
    Aggregated counts of tests per experiment
    :param pids: (list) of files to be aggregated
    :param output_file: (str) of output file location
    :param cl: (boolean) if cell_line should be considered or not
    :return: -
    """
    ppi_dict = defaultdict(lambda: [0, 0, set()])
    bait_idx = dict()
    prey_idx = dict()
    cl_idx = dict()
    pid_idx = dict()

    for i, study in enumerate(pids):
        with (open(study, "r") as f):
            next(f)  # header
            for line in f:
                values = line.strip().split("\t")
                bait, prey, n_tested, n_observed = values[:4]
                pid = values[-2]
                keys = (update_idx(bait_idx, bait), update_idx(prey_idx, prey), update_idx(cl_idx, values[-1]))
                pid_id = update_idx(pid_idx, pid)
                ppi_dict[keys][0] += int(n_tested)
                ppi_dict[keys][1] += int(n_observed)
                ppi_dict[keys][2] |= {pid_id}

    bait_idx = {value: key for key, value in bait_idx.items()}
    prey_idx = {value: key for key, value in prey_idx.items()}
    cl_idx = {value: key for key, value in cl_idx.items()}
    pid_idx = {value: key for key, value in pid_idx.items()}


    with open(output_file, "w") as w:
        header_line = "gene_name_bait\tgene_name_prey\tn_tested\tn_observed\tpubmed_id\tcl_id\n"
        w.write(header_line)
        for keys, observed_tested_pids in ppi_dict.items():
            pids = ";".join(sorted(
                [pid_idx[pid] for pid in observed_tested_pids[2]]))
            line = "\t".join([bait_idx[keys[0]], prey_idx[keys[1]]]) + "\t"
            line += "\t".join(map(str, observed_tested_pids[:2])) + "\t"
            line += pids + "\t"
            line += cl_idx[keys[2]] + "\n"
            w.write(line)

