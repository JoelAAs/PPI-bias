configfile: "config_files/config.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments


#### Config
pn = config["project_name"]
if pn:
    pn = "/" + pn

datasets = [
    "flat",
    "ms",
    "y2h",
    #"abundance_mcmc",
    #"MI-1314"
]

pods = [
    f"work_folder{pn}/analysis/POD/POD_{data}.csv" for data in datasets
]

colocalisation_plot = [
    f"work_folder{pn}/plots/AccumulationPOD/colocalisation_{data}.png"
    for data in datasets
]

matched_colocalisation_plot = [
    f"work_folder{pn}/plots/AccumulationPOD/matched_colocalisation_{data}.png"
    for data in datasets if data != "abundance_mcmc"
]

go_jaccards_plot = [
    f"work_folder{pn}/plots/AccumulationPOD/go_{data}_jaccard.png"
    for data in datasets
]

do_jaccards_plot = [
    f"work_folder{pn}/plots/AccumulationPOD/do_{data}_jaccard.png"
    for data in datasets
]

hydro_delta_plot = [
    f"work_folder{pn}/plots/AccumulationPOD/hydrophobicity_{data}.png"
    for data in datasets
]

### Negatome compare
negatome_compare = [
    f"work_folder{pn}/analysis/neg2compare/{data}.txt"
    for data in datasets
]
negatome_entropy = [
    f"work_folder{pn}/analysis/negatome/test_entropy_{data}_limit_{min_tests}.csv"
    for min_tests in [3, 4, 5] for data in ["y2h", "ms"]
]

expected_output = pods
#expected_output += colocalisation_plot + go_jaccards_plot + hydro_delta_plot +  do_jaccards_plot
#expected_output += negatome_compare + matched_colocalisation_plot + negatome_entropy


## Sub workflows
include: "src/FormatFiltering/FormatingFiltering.smk"

include: "src/ExperimentalSearchSpace/experimental_search_space.smk"
include: "src/ExperimentalSearchSpace/CountProteinPairs.smk"

include: "src/Analysis/CellLine/cell_line_analysis.smk"
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
    model_configuration="[a-z0-9]+"

rule all:
    input:
        expand(
            f"work_folder{pn}/subsets/report/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}.nb.html",
            dataset=datasets,neg_limit=[2,3],pos_limit=0.15,partition_name=["sequencesimilarity", "maxpos"]
        ),
        expand(
            f"work_folder{pn}/classification/randomforest/{{dataset}}_{{model_configuration}}_model_parameters.txt",
            dataset=["ms", "y2h", "flat"], model_configuration = ["seqs2","seqs3", "maxdata2", "maxdata3"]),
        f"work_folder{pn}/classification/randomforest/goldensplit_asis_model_parameters.txt"

    #expected_output,,
    #f"work_folder{pn}/embeddings/canonical_embedding.csv.gz",
    #f"work_folder{pn}/plots/degree/GO_enrichment.png",
    #f"work_folder{pn}/plots/localisation/HuRI_bioplex.png",
    #f"work_folder{pn}/plots/membrane/HuRI_bioplex.png"
