configfile: "config_files/config_uniport_pod.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments
from src.support_functions import read_fasta

#### Config

datasets = config["datasets"]
esm_models = ["ESM2", "ESMC"]
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

include: "src/Analysis/Annotation/AnnotationProbabilities.smk"

include: "src/Analysis/BiologicalInsights/Shared.smk"
include: "src/Analysis/Expression/Coexpression.smk"
include: "src/Analysis/BiologicalInsights/AssayConcordance.smk"
include: "src/Analysis/BiologicalInsights/DegreeQuadrants.smk"
include: "src/Analysis/BiologicalInsights/NoInteractionHubs.smk"

include: "src/Analysis/NegatomeComparison/NegatomeAnalysis.smk"
include: "src/Analysis/NegatomeComparison/CompareSharedBaits.smk"
include: "src/Analysis/CompareLocalisationMethod/MethodLocalisation.smk"

include: "src/PPIClassification/Embeddings/Embeddings.smk"
include: "src/PPIClassification/DataSplit/GetGraphs.smk"
include: "src/PPIClassification/DataSplit/GenePartitions.smk"
include: "src/PPIClassification/DataSplit/GenerateSplits.smk"
include: "src/PPIClassification/DataSplit/GetGoldenSplit.smk"
include: "src/PPIClassification/DataSplit/CheckRedundancy.smk"
include: "src/PPIClassification/ModelEvaluation/Evaluations.smk"


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
    classifier="[a-z]+"

rule all:
    input:
        expand("work_folder/classification/{classifier}/permuted/all_metrics_{network_type}_{esm_model}.csv",
            classifier="xgboost", network_type="undirectional",esm_model="ESM2"),
        #"work_folder/analysis/shared_annotation_proportions/plots/undirectional_OR.png",
        # Module 1: co-expression
        "work_folder/analysis/coexpression/plots/undirectional_coexpression.png",
        "work_folder/analysis/coexpression/plots/undirectional_or_high_membership.png",
        "work_folder/analysis/coexpression/plots/undirectional_summed_coexpression_vs_random.png",
        # Module 4: assay concordance
        "work_folder/analysis/assay_concordance/plots/concordance.png",
        "work_folder/analysis/protein_degree/plots/undirectional_degree_distributions.png",
        "work_folder/analysis/protein_degree/plots/undirectional_degree_bin_heatmap.png",
        # Module 3: degree quadrants
        "work_folder/analysis/degree_quadrants/plots/undirectional_quadrants.png",
        "work_folder/analysis/degree_quadrants/plots/undirectional_degree_diff_vs_references.png",
        #expand("work_folder/analysis/degree_quadrants/{dataset}_undirectional_go_enrichment.tsv",
        #    dataset=config["datasets"]),
        ## Module 2: no-interaction hubs
        "work_folder/analysis/no_interaction_hubs/plots/undirectional_hub_properties.png"
