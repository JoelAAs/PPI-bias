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
include: "src/Analysis/Localisation/HCI_Negatome_test.smk"
include: "src/Plotting/get_plots.smk"
include: "src/Analysis/AbundanceAwareDetection/MCMC_abundance.smk"
include: "src/Analysis/Enrichment/GetDegree.smk"
include: "src/Analysis/Enrichment/EnrichmentGODO.smk"
include: "src/Analysis/PairFunctionality/GetSharedFunctionality.smk"

wildcard_constraints:
    cell_line="_[_a-zA-Z]+"

# expected_output = [
#     f"work_folder/inferred_search_space/aggregated/multi_methods/{multi_method}_experimental_wise.csv" for
#         multi_method in ["ms", "y2h"]
# ]

expected_output = [
    "work_folder/analysis/abundance_aware/bait_prey_abundance.csv",
    "work_folder/analysis/abundance_aware/localisation/probability_match_abundance_less.csv",
    "work_folder/inferred_search_space/analysis/bias_reduced_ppis/localisation_p_estimated_protein_pairs_less.csv"
    "work_folder/plots/degree/doid_vs_deg.png",
    "work_folder/analysis/GO/flat_jaccard.csv",
    "work_folder/analysis/GO/abundance_jaccard.csv"

]

rule all:
    input:
        expected_output


