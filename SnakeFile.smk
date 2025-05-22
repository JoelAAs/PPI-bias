#configfile: "config_cell.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments


include: "src/FormatFiltering/FormatingFiltering.smk"
include: "src/ExperimentalSearchSpace/experimental_search_space.smk"
include: "src/Analysis/CellLine/cell_line_analysis.smk"
include: "src/Analysis/DetectionMethod/detection_method.smk"
include: "src/Analysis/Localisation/localisation.smk"
include: "src/Analysis/ExperimentalNegatome/experimental_negatome.smk"
include: "src/Plotting/get_plots.smk"


expected_output = [
    f"work_folder/inferred_search_space/aggregated/multi_methods/{multi_method}_experimental_wise.csv" for
        multi_method in ["ms", "y2h"]
]

expected_output += [
    "work_folder/plots/localisation_OR_y2h_ms.png",
    "work_folder/inferred_search_space/aggregated/methods/threshold_negatome.csv",
    "work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"
    "work_folder/plots/cell_line_prey.png"
    "work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
]

rule all:
    input:
        expected_output