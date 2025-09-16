#configfile: "config_cell.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments


include: "src/FormatFiltering/FormatingFiltering.smk"
include: "src/ExperimentalSearchSpace/experimental_search_space.smk"
include: "src/Analysis/CellLine/cell_line_analysis.smk"
include: "src/Analysis/DetectionMethod/detection_method.smk"
include: "src/Analysis/Annotation/CoLocalisation.smk"
include: "src/Analysis/ExperimentalNegatome/experimental_negatome.smk"
include: "src/Plotting/get_plots.smk"
include: "src/Analysis/AbundanceAwareDetection/MCMC_abundance.smk"
include: "src/Analysis/Enrichment/GetDegree.smk"
include: "src/Analysis/Enrichment/EnrichmentGODO.smk"
include: "src/Analysis/Annotation/OverlapGO.smk"
include: "src/Analysis/Annotation/HydrophobicitySimilarity.smk"

wildcard_constraints:
    cell_line="_[_a-zA-Z]+"

datasets = [
    "ms_y2h",
    "ms",
    "y2h",
    "abundance",
    "MI-1314"
]
pods = [
    f"work_folder/analysis/POD/POD_{data}.csv" for data in datasets
]

expected_output = pods

rule all:
    input:
        expected_output


