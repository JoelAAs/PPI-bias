configfile: "config_test.yaml"

include: "src/FormatingFiltering.smk"
include: "src/ProbabilityOfPreyDetection/ProbabilityOfPreyDetection.smk"
include: "src/CellLineFunctions/GetCellLineFromPID.smk"
include: "src/CellLineFunctions/CellLinePPIOverlap.smk"
include: "src/CellLineFunctions/GetPubtatorCount.smk"

method_aggregate = dict()
for method in config["methods"]:
    try:
        config[config["methods"][method]].append(method)
    except KeyError:
        config[config["methods"][method]] = [method,]




rule all:
    input:
        # "work_folder/intact/pair_count/ppi_pair_counts_all_MI-1112.csv",
        # "work_folder/ppi_cl_overlap/one_hot_MI-0007.csv"
        "work_folder/ppi_cl_overlap/IoU_MS.csv",
        "work_folder/intact/pubtator_cl_count/overlap/MI-0007.csv"
        # expand(
        #     "work_folder/intact/pair_count/ppi_pair_counts_{cell_line}_{method}.csv",
        #     method = config["methods"], cell_line = "HeLa"
        # )
