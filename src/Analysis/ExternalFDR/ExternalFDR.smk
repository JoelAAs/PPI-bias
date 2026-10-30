import functools
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc


@functools.cache
def _pod_study_ids(network_type, dataset):
    pod = f"work_folder/analysis/POD/{network_type}/POD_{dataset}.pq"
    prefix = workflow.storage_settings.default_storage_prefix
    pod_path = f"{prefix.rstrip('/')}/{pod}" if prefix else pod
    col = pq.read_table(pod_path, columns=["pubmed_id"], read_dictionary=["pubmed_id"])["pubmed_id"]
    joined = pc.unique(pa.chunked_array([c.dictionary for c in col.chunks]))
    return sorted(pc.unique(pc.list_flatten(pc.split_pattern(joined, ";"))).to_pylist())


def get_consituent_studies(wc):
    """
    The studies POD was aggregated from, read back out of its pubmed_id column
    :param wc: snakemake wildcards
    :return: (list) one inferred search space file per study
    """
    study_folder = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]
    study_folder = "work_folder" + study_folder.split("work_folder")[1]

    studies = _pod_study_ids(wc.network_type, wc.dataset)
    return expand(f"{study_folder}/{{study}}.csv", study=studies)




rule leave_one_out_FDR_TPR:
    """
    Leave each study out of POD in turn and score the pairs it reported against the rest
    """
    params:
        id_pattern = config["id_pattern"],
        n_em_iterations = config["n_em_iterations"],
        shrinkage_pseudocount = config["shrinkage_pseudocount"],
        n_reference_studies = config["n_reference_studies"],
    input:
        pod_file = "work_folder/analysis/POD/directional/POD_{dataset}.pq",
        consituent_studies = get_consituent_studies
    output:
        study_fdr = "work_folder/analysis/FDR_aware/study_metrics/directional_{dataset}.tsv",
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/directional_{dataset}.pq",
        global_fit = "work_folder/analysis/FDR_aware/global_fit/directional_{dataset}.tsv",
    threads: 20
    log:
        "logs/analysis/FDR_aware/study_metrics/directional_{dataset}.log"
    script:
        "scripts/per_pair_study_fdr.py"

rule leave_one_out_FDR_TPR_directional:
    input:
        study_fdr = "work_folder/analysis/FDR_aware/study_metrics/directional_{dataset}.tsv",
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/directional_{dataset}.pq",
        global_fit = "work_folder/analysis/FDR_aware/global_fit/directional_{dataset}.tsv",
    output:
        study_fdr = "work_folder/analysis/FDR_aware/study_metrics/undirectional_{dataset}.tsv"
    # TODO:!


rule plot_log_odds_vs_pod:
    input:
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/{network_type}_{dataset}.pq",
        pod_file = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq"
    output:
        plot = "work_folder/analysis/FDR_aware/plots/{network_type}_{dataset}_log_odds_vs_pod.png",
        plot_not_observed = "work_folder/analysis/FDR_aware/plots/{network_type}_{dataset}_log_odds_vs_pod_not_observed.png"
    log:
        "logs/analysis/FDR_aware/plots/{network_type}_{dataset}_log_odds_vs_pod.log"
    script:
        "scripts/plot_log_odds_vs_pod.R"


rule plot_observable_prey_vs_bait_count:
    params:
        ms_methods = config["ms"]
    input:
        bp_df = "work_folder/formated/bait_prey_publications.csv"
    output:
        plot = "work_folder/analysis/FDR_aware/plots/bait_vs_prey.png",
        fit_stats = "work_folder/analysis/FDR_aware/plots/bait_vs_prey_fit.tsv"
    script:
        "scripts/baits_vs_bait_saturation.R"


