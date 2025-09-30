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
include: "src/Analysis/Annotation/OverlapDO.smk"
include: "src/Analysis/Annotation/HydrophobicitySimilarity.smk"
include: "src/Analysis/NegatomeComparison/NegatomeAnalysis.smk"

wildcard_constraints:
    cell_line="_[_a-zA-Z]+",
    subset="[a-zA-Z0-9-]+",
    model="[_a-zA-Z0-9-]+",
    pid="[:a-zA-Z0-9-]+"

datasets = [
    "flat",
    "ms",
    "y2h",
    "abundance_mcmc",
    "MI-1314"
]
pods = [
    f"work_folder/analysis/POD/POD_{data}.csv" for data in datasets
]

### annotations
colocalisation = [
    f"work_folder/analysis/localisation/cumulative/POD_{data}_localisation_lesser.csv"
    for data in datasets
]
go_jaccards = [
    f"work_folder/analysis/GO/cumulative/POD_{data}_jaccard_lesser.csv"
    for data in datasets
]

hydro_delta = [
    f"work_folder/analysis/hydrophobicity/cumulative/POD_{data}_netsurfp2_lesser.csv"
    for data in datasets
]
### Plotting
colocalisation_plot = [
    f"work_folder/plots/AccumulationPOD/colocalisation_{data}.png"
    for data in datasets
]

matched_colocalisation_plot = [
    f"work_folder/plots/AccumulationPOD/matched_colocalisation_{data}.png"
    for data in datasets
]

go_jaccards_plot = [
    f"work_folder/plots/AccumulationPOD/go_{data}_jaccard.png"
    for data in datasets
]

do_jaccards_plot = [
    f"work_folder/plots/AccumulationPOD/do_{data}_jaccard.png"
    for data in datasets
]

hydro_delta_plot = [
    f"work_folder/plots/AccumulationPOD/hydrophobicity_{data}.png"
    for data in datasets
]

### Negatome compare
negatome_compare = [
    f"work_folder/analysis/neg2compare/{data}.txt"
    for data in datasets
]


expected_output = pods + colocalisation + go_jaccards + hydro_delta
expected_output += colocalisation_plot + go_jaccards_plot + hydro_delta_plot +  do_jaccards_plot
expected_output += negatome_compare + matched_colocalisation_plot

rule all:
    input:
        #expected_output,
        matched_colocalisation_plot,
        "work_folder/plots/degree/GO_enrichment.png"


