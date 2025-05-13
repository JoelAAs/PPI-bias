import pandas as pd

# Run only with the correct config file
# It won't run otherwise

from collections import defaultdict

def nested_dict():
    return defaultdict(nested_dict)

def get_input_for_aggregation(wc , filename):
    CL_FOLDER = checkpoints.infer_experimental_search_space.get().output[0]
    cl_df = pd.read_csv(filename, sep="\t")
    cl_df = cl_df[["pubmed_id", "detection_method", "cl_id"]]
    cl_df = cl_df[~cl_df.duplicated()]

    expected_input = [
        f"{CL_FOLDER}/{pubmed_id}_{detection_method}_{cl_id}.csv"
        for _, (pubmed_id, detection_method, cl_id) in cl_df.iterrows()
    ]
    return expected_input


rule aggregate_inferred_studies_cell_line:
    input:
        cl_pids = lambda wc: get_input_for_aggregation(wc, config["ppi_df"])
    output:
        cell_line_counts = "work_folder/inferred_search_space/aggregated/cell_line_specific.csv"
    run:
        ppi_dict = nested_dict()
        for cl_study in input.cl_pids:
            with open(cl_study, "r") as f:
                header = True
                for line in f:
                    if header:
                        header = False
                    else:
                        bait, prey ,n_tested, n_observed, detection_method, pubmed_id, cl_id = line.strip().split("\t")

                        if not ppi_dict[bait][prey][cl_id]:
                            ppi_dict[bait][prey][cl_id]["n_tested"] = 0
                            ppi_dict[bait][prey][cl_id]["n_observed"] = 0

                        ppi_dict[bait][prey][cl_id]["n_tested"] += int(n_tested)
                        ppi_dict[bait][prey][cl_id]["n_observed"] += int(n_observed)

        with open(output.cell_line_counts, "w") as w:
            w.write("gene_name_bait\tgene_name_prey\tn_tested\tn_observed\tcl_id\n")

            for c_bait, prey_dict in ppi_dict.items():
                for c_prey, cl_dict in prey_dict.items():
                    for c_cl, count_dict in cl_dict.items():
                        w.write(
                            "\t".join([
                                c_bait,
                                c_prey,
                                str(count_dict["n_tested"]),
                                str(count_dict["n_observed"]),
                                c_cl
                            ]) + "\n"
                        )



