def get_consituent_studies(wc):
    """
    The studies POD was aggregated from, read back out of its pubmed_id column
    :param wc: snakemake wildcards
    :return: (list) one inferred search space file per study
    """
    study_folder = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]
    study_folder = "work_folder" + study_folder.split("work_folder")[1]

    pod = f"work_folder/analysis/POD/{wc.network_type}/POD_{wc.dataset}.pq"
    prefix = workflow.storage_settings.default_storage_prefix
    storage_file = storage.fs(f"{prefix.rstrip('/')}/{pod}" if prefix else pod)
    storage_object = storage_file.flags["storage_object"]
    storage_object.local_path().parent.mkdir(parents=True, exist_ok=True)
    storage_object.retrieve_object()

    pids = pd.read_parquet(storage_file, columns=["pubmed_id"])["pubmed_id"]
    studies = {p for ids in pids for p in ids.split(";")}
    return expand(f"{study_folder}/{{study}}.csv", study=sorted(studies))



rule global_detection_mixture:
    params:
        n_em_iterations = 10000
    input:
        pod_file = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
    output:
        study_fdr = "work_folder/analysis/FDR_aware/study_metrics/global_{network_type}_{dataset}.tsv"
    script:
        "scripts/estimate_global_fdr.py"


rule leave_one_out_FDR_TPR:
    """
    Leave each study out of POD in turn and score the pairs it reported against the rest
    """
    params:
        id_pattern = config["id_pattern"],
        n_em_iterations = config["n_em_iterations"]
    input:
        pod_file = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
        consituent_studies = get_consituent_studies
    output:
        study_fdr = "work_folder/analysis/FDR_aware/study_metrics/{network_type}_{dataset}.tsv",
        pair_posterior = "work_folder/analysis/FDR_aware/pair_posterior/{network_type}_{dataset}.pq"
    threads: 20
    log:
        "logs/analysis/FDR_aware/study_metrics/{network_type}_{dataset}.log"
    script:
        "scripts/per_pair_study_fdr.py"


rule plot_observable_prey_vs_bait_count:
    params:
        ms_methods = config["ms"]
    input:
        bp_df = "work_folder/formated/bait_prey_publications.csv"
    output:
        plot = "work_folder/analysis/FDR_aware/plots/bait_vs_prey.png"
    script:
        "scripts/baits_vs_bait_saturation.R"