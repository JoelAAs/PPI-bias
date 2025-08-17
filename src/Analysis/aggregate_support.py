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
    ppi_dict = defaultdict(lambda: [0, 0])
    bait_idx = prey_idx = cl_idx = dict()
    for i, study in enumerate(pids):
        with (open(study, "r") as f):
            next(f)  # header
            for line in f:
                values = line.strip().split("\t")
                bait, prey, n_tested, n_observed = values[:4]
                keys = (update_idx(bait_idx, bait), update_idx(prey_idx, prey)) if not cl else (
                    update_idx(bait_idx, bait), update_idx(prey_idx, prey), update_idx(cl_idx, values[-1]))
                ppi_dict[keys][0] += int(n_tested)
                ppi_dict[keys][1] += int(n_observed)

    bait_idx ={value: key for key, value in bait_idx.items()}
    prey_idx = {value: key for key, value in prey_idx.items()}
    cl_idx = {value: key for key, value in cl_idx.items()}

    with open(output_file, "w") as w:
        if cl:
            header_line = "gene_name_bait\tgene_name_prey\tcl_id\tn_tested\tn_observed\n"
        else:
            header_line = "gene_name_bait\tgene_name_prey\tn_tested\tn_observed\n"
        w.write(header_line)
        for keys, observed_tested in ppi_dict.items():
            protein_pair = "\t".join([bait_idx[keys[0]], prey_idx[keys[1]]])
            if cl:
                protein_pair += "\t" + cl_idx[keys[2]]
            w.write(protein_pair + "\t" + "\t".join(map(str, observed_tested)) + "\n")
