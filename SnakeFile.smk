configfile: "config_files/config.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments
from src.support_functions import read_fasta

#### Config
pn = config["project_name"]
if pn:
    pn = "/" + pn

datasets = config["datasets"]
#expected_output += colocalisation_plot + go_jaccards_plot + hydro_delta_plot +  do_jaccards_plot
#expected_output += negatome_compare + matched_colocalisation_plot + negatome_entropy


## Sub workflows
include: "src/FormatFiltering/FormatingFiltering.smk"

include: "src/ExperimentalSearchSpace/experimental_search_space.smk"
include: "src/ExperimentalSearchSpace/CountProteinPairs.smk"

include: "src/Analysis/CellLine/cell_line_analysis.smk"
include: "src/Analysis/CellLine/GetCellLineAnnotation.smk"
include: "src/Analysis/DetectionMethod/detection_method.smk"
include: "src/Analysis/ExperimentalNegatome/experimental_negatome.smk"
include: "src/Analysis/AbundanceAwareDetection/MCMC_abundance.smk"

include: "src/Analysis/Enrichment/GetDegree.smk"
include: "src/Analysis/Enrichment/EnrichmentAnalysisGeneSet.smk"

include: "src/Analysis/Annotation/CoLocalisation.smk"
include: "src/Analysis/Annotation/OverlapGO.smk"
include: "src/Analysis/Annotation/OverlapDO.smk"
include: "src/Analysis/Annotation/HydrophobicitySimilarity.smk"
include: "src/Analysis/Annotation/InterfaceStatistics.smk"

include: "src/Analysis/NegatomeComparison/NegatomeAnalysis.smk"
include: "src/Analysis/NegatomeComparison/CompareSharedBaits.smk"
include: "src/Analysis/CompareLocalisationMethod/MethodLocalisation.smk"

include: "src/PPIClassification/Embeddings/Embeddings.smk"
include: "src/PPIClassification/DataSplit/GetGraphs.smk"
include: "src/PPIClassification/DataSplit/GenePartitions.smk"
include: "src/PPIClassification/DataSplit/GenerateTrainTestSplits.smk"
include: "src/PPIClassification/DataSplit/BalanceSplits.smk"
include: "src/PPIClassification/DataSplit/GetGoldenSplit.smk"
include: "src/PPIClassification/DataSplit/CheckRedundancy.smk"
include: "src/PPIClassification/ModelEvaluation/Evaluations.smk"


include: "src/PPIClassification/Classification/RandomForest.smk"
include: "src/PPIClassification/Report/Reporting.smk"
include: "src/Plotting/get_plots.smk"

wildcard_constraints:
    cell_line="_[_a-zA-Z]+",
    subset="[a-zA-Z0-9-]+",
    model="[_a-zA-Z0-9-]+",
    data="[_a-zA-Z0-9-]+",
    dataset="[_a-zA-Z0-9-]+",
    pid="[:a-zA-Z0-9-]+",
    neg_limit="[0-9.]+",
    pos_limit="[0-9.]+",
    model_configuration="[a-z0-9]+",
    selected_data="[a-z0-9_.]+",
    network_type="(directional|undirectional)"

rule all:
    input:
        expand("work_folder{pn}/subsets/train/equal_edge/{dataset}_directional_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
            pn=pn, dataset=datasets, pos_limit=config["positive_limits"], neg_limit=config["negative_limits"])
