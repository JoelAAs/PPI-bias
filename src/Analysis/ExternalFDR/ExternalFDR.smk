import functools
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc


@functools.cache
def _pod_study_ids(dataset):
    pod = f"work_folder/analysis/POD/directional/POD_{dataset}.pq"
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
    pod_path = checkpoints.all_methods_filter_out.get(network_type="directional", data=wc.dataset).output.full_detection # need to exist first

    studies = _pod_study_ids(wc.dataset)
    return expand(f"{study_folder}/{{study}}.csv", study=studies)




rule leave_one_out_FDR_TPR_directional:
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
        refits = "work_folder/analysis/FDR_aware/diagnostics/refits_directional_{dataset}.tsv",
        predictive_check = "work_folder/analysis/FDR_aware/diagnostics/predictive_check_directional_{dataset}.tsv",
        floor_sensitivity = "work_folder/analysis/FDR_aware/diagnostics/floor_sensitivity_directional_{dataset}.tsv",
    threads: 20
    log:
        "logs/analysis/FDR_aware/study_metrics/directional_{dataset}.log"
    script:
        "scripts/per_pair_study_fdr.py"

rule leave_one_out_FDR_TPR_undirectional:
    params:
        id_pattern = config["id_pattern"],
    input:
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/directional_{dataset}.pq",
        global_fit = "work_folder/analysis/FDR_aware/global_fit/directional_{dataset}.tsv",
    output:
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/undirectional_{dataset}.pq"
    script:
        "scripts/join_directional.py"


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


