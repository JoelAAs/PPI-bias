configfile: "config_test.yaml"

include: "src/FormatingFiltering.smk"
include: "src/ProbabilityOfPreyDetection/ProbabilityOfPreyDetection.smk"
include: "src/CellLineExtraction/GetCellLineFromPID.smk"

rule all:
    input:
        expand(
            "work_folder/intact/pair_count/ppi_pair_counts_{cell_line}_{method}.csv",
            method = config["methods"], cell_line = "HeLa"
        )
