configfile: "config_files/config_uniport_pod.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments
from src.support_functions import read_fasta

#### Config

datasets = config["datasets"]
esm_models = ["ESM2", "ESMC"]
prediction_models = ["AF", "RF2-PPI", "DCA"]
network_types = ["directional", "undirectional"]
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
include: "src/Analysis/ExternalFDR/ExternalFDR.smk"
include: "src/Analysis/AbundanceAwareDetection/MCMC_abundance.smk"

include: "src/Analysis/Annotation/AnnotationProbabilities.smk"
include: "src/Analysis/Annotation/InterfaceStatistics.smk"
include: "src/Analysis/InterfaceStatistics/ProteinInterfaces.smk"

include: "src/Analysis/Expression/Coexpression.smk"
include: "src/Analysis/Co-essentiality/CoEssensiality.smk"
include: "src/Analysis/NonInteractionHubs/NoInteractionHubs.smk"
include: "src/Analysis/NonInteractionHubs/DegreeQuadrants.smk"
include: "src/Analysis/MethodConcordance/MethodConcordance.smk"
include: "src/Analysis/ProteinPromiscuity/ProteinPromiscuity.smk"
include: "src/Analysis/ProteinPromiscuity/AutoActivatorsContaminants.smk"

include: "src/Analysis/NegatomeComparison/NegatomeAnalysis.smk"
include: "src/Analysis/NegatomeComparison/CompareSharedBaits.smk"
include: "src/Analysis/CompareLocalisationMethod/MethodLocalisation.smk"
include: "src/Analysis/CheckAgainstPredictions/CheckAgainstPredicted.smk"
include: "src/Analysis/CheckAgainstPredictions/CheckAgainstOtherMethods.smk"

include: "src/PPIClassification/Embeddings/Embeddings.smk"
include: "src/PPIClassification/DataSplit/GetGraphs.smk"
include: "src/PPIClassification/DataSplit/GenePartitions.smk"
include: "src/PPIClassification/DataSplit/GenerateSplits.smk"
include: "src/PPIClassification/DataSplit/GetGoldenSplit.smk"
include: "src/PPIClassification/DataSplit/CheckRedundancy.smk"
include: "src/PPIClassification/ModelEvaluation/Evaluations.smk"
include: "src/PPIClassification/DatasetBalance/DatasetBalance.smk"


include: "src/PPIClassification/Classification/Classifiers.smk"
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
    pos_limit="([0-9.]+|all)",
    model_configuration="[a-z0-9]+",
    selected_data="[a-z0-9_.]+",
    network_type="(directional|undirectional)",
    random="(-random)?",
    esm_model="[A-Z0-9]+",
    permutation="[0-9]+",
    classifier="[a-z]+",
    pair_set="(all|correct)"

rule all:
    input:
        expand("work_folder/analysis/POD/{network_type}/POD_{dataset}.pq", dataset=datasets, network_type = network_types),
        expand("work_folder/analysis/FDR_aware/study_metrics/{network_type}_{dataset}.tsv", dataset=datasets, network_type = network_types)