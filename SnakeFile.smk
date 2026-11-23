configfile: "config_files/config_uniport_pod.yaml"
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact, false_discovery_control
from src.Analysis.aggregate_support import aggregate_inferred_experiments
from src.support_functions import read_fasta

#### Config

datasets = config["datasets"]
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
include: "src/PPIClassification/DataSplit/NaiveSplit.smk"
include: "src/PPIClassification/DataSplit/CheckRedundancy.smk"
include: "src/PPIClassification/ModelEvaluation/Evaluations.smk"
include: "src/PPIClassification/DatasetBalance/DatasetBalance.smk"


include: "src/PPIClassification/Classification/Classifiers.smk"
include: "src/PPIClassification/Classification/EvaluatePredefinedSplits.smk"
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
    esm_model="[A-Za-z0-9]+",
    permutation="[0-9]+",
    classifier="[a-z]+",
    pair_set="(all|correct)",
    source="[a-zA-Z0-9]+",
    pos_set="(hippie|fdr|pod)",
    neg_set="(fdr|pod)",
    no="(-no)?",
    test_no="(-no)?",
    seed="[0-9]+",
    edge_type="(pos|neg)",
    exp_set="(fdr|pod)",
    set="seed_[0-9]+_(train|test)",
    state="(cand|nonobs)",
    test_state="(cand|nonobs)"

rule all:
    input:
        # expand("work_folder/classification/{source}/{classifier}/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.png",
        #        **config["classification_targets"]),
        # expand("work_folder/classification/{source}/{classifier}/hrni_no_delta/{network_type}_{esm_model}_delta_auc_summary.tsv",
        #        **config["classification_targets"]),
        # expand("work_folder/analysis/FDR_aware/global_fit/{network_type}_{dataset}.tsv", dataset=["ms"], network_type = "directional"),
        # "work_folder/analysis/FDR_aware/plots/undirectional_ms_fdr_vs_pod.png",
        # expand("work_folder/analysis/ProteinPromiscuity/autoactivators/groups/{dataset}_feature_log_or.svg",
        #     dataset=["ms", "y2h"]),
        # expand("work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{dataset}_detection_distribution.png", dataset=["ms", "y2h"]),  
        # expand("work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{network_type}_{dataset}_pod_bait_prey_degree_correlation.png",
        #     network_type=network_types, dataset=["ms", "y2h"]),
        # "work_folder/classification/plots/hrni_no_auc_arrows.png",
        #expand("work_folder/classification/naive_split/{classifier}/all_metrics_{network_type}_{esm_model}.tsv",
        #    classifier=config["naive_split"]["classifier"], network_type=config["naive_split"]["network_type"],
        #    esm_model=config["naive_split"]["esm_model"]),
        #"work_folder/classification/predef/xgboost/all_metrics.tsv",
        expand("work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz", esm_model=["ESM2", "ProstT5"]),

        #expand("work_folder/classification/naive_split/plots/{network_type}_{classifier}_{plot}.svg",
        #    network_type=config["naive_split"]["network_type"], classifier=config["naive_split"]["classifier"],
        #    plot=["auc_hrni_test", "auc_hrni_test_by_embedding", "delta_hrni_vs_no_test", "delta_auc_vs_adjusted"]),
        #"work_folder/classification/embedding_sim_auc/all_classifications.tsv",
        # expand("work_folder/classification/embedding_sim_auc/plots/{dataset}_{source}_{network_type}_{esm_model}_cosine_sim_vs_score.png",
        #     source = ["pod", "fdr"], dataset="ms", network_type=["directional", "undirectional"], esm_model=["ESM2", "ProstT5"]
        # ),
        # expand("work_folder/classification/embedding_sim_auc/plots/{dataset}_{source}_{network_type}_{esm_model}_cosine_sim_vs_score.png",
        #     source = "pod", dataset="y2h", network_type=["directional", "undirectional"], esm_model=["ESM2", "ProstT5"]
        # ),
        #expand("work_folder/classification/{classifier}/hrni_no_delta/plots/{network_type}_{esm_model}_delta_auc_boxplot.png",
        #       classifier="xgboost", network_type=network_types, esm_model=["ESM2", "ProstT5"]),
        #expand("work_folder/analysis/FDR_aware/plots/{network_type}_ms_log_odds_vs_pod.png", network_type=network_types),
        #"work_folder/analysis/ProteinPromiscuity/autoactivators/known/ms_prey_added.pdf",
        #"work_folder/analysis/ProteinPromiscuity/autoactivators/known/y2h_bait_added.pdf",
        #expand("work_folder/analysis/ProteinPromiscuity/features/tests/{dataset}_{network_type}_{feature}.tsv",
        # dataset=["ms", "y2h"], network_type=["directional", "undirectional"], feature = ["nucleoplasm-vs-cytocol",]),