configfile: "config_test.yaml"

include: "src/FormatingFiltering.smk"
include: "src/ProbabilityOfPreyDetection/ProbabilityOfPreyDetection.smk"
include: "src/CellLineFunctions/GetCellLineFromPID.smk"
include: "src/CellLineFunctions/CellLinePPIOverlap.smk"

current_ms_cl =[
    "HeLa",
    "HEK293",
    "HEK293T",
    "U2OS",
]

rule all:
    input:
        expand(
            "work_folder/intact/pair_count/ppi_pair_counts_{cl}_MI-0007.csv",
            cl = current_ms_cl),
        "work_folder/intact/pair_count/ppi_pair_counts_all_MI-1112.csv",
        "work_folder/ppi_cl_overlap/one_hot_MI-0007.csv"
        # expand(
        #     "work_folder/intact/pair_count/ppi_pair_counts_{cell_line}_{method}.csv",
        #     method = config["methods"], cell_line = "HeLa"
        # )
