#configfile: "config_cell.yaml"

config["ppi_df"] = "../Cl-annotated-ppis/work_folder/CL_annotated_bait_prey.csv"
config["id_pattern"] = "gene_name"
config["cell_line_present"] = True

include: "src/FormatFiltering/FormatingFiltering.smk"
include: "src/ExperimentalSearchSpace/experimental_searchspace.smk"
include: "src/Analysis/CellLineAnalysis/cell_line_analysis.smk"


rule all:
    input:
        "work_folder/inferred_search_space/aggregated/cell_line_specific.csv"
