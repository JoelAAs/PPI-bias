import random

def select_row_wise_study_files(study_folder, methods, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = study_folder

    # storage.fs() is not auto-prefixed, so add the prefix and retrieve the file ourselves
    prefix = workflow.storage_settings.default_storage_prefix
    query = f"{prefix.rstrip('/')}/{filename}" if prefix else filename
    storage_file = storage.fs(query)
    storage_object = storage_file.flags["storage_object"]
    storage_object.local_path().parent.mkdir(parents=True, exist_ok=True)
    storage_object.retrieve_object()
    ppi_df = pd.read_csv(storage_file, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"].isin(methods)]
    ppi_df = ppi_df[
        ppi_df[f"{id_pattern}_bait"] != ppi_df[f"{id_pattern}_prey"]
    ]

    if ppi_df.empty:
        raise ValueError(f"No studies for method: {methods}")
    ppi_df = ppi_df[~ppi_df[[f"{id_pattern}_bait", f"{id_pattern}_prey", "pubmed_id"]].duplicated(keep="first")] # Remove isoforms, if considering gene names
    if remove_single:
        counts = ppi_df.groupby(["pubmed_id"], as_index=False).size()
        keep_pids = counts[counts["size"] != 1]["pubmed_id"]
        ppi_df = ppi_df[ppi_df["pubmed_id"].isin(keep_pids)]
    expected = [
        f"{STUDY_FOLDER}/{pid}_{method}.csv"
        for pid, method in ppi_df[["pubmed_id", "detection_method"]].drop_duplicates().itertuples(index=False)
    ]
    return expected

def get_row_wise_input_files(methods, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]
    STUDY_FOLDER = "work_folder" + STUDY_FOLDER.split("work_folder")[1]
    return select_row_wise_study_files(STUDY_FOLDER, methods, id_pattern, filename, remove_single)

rule get_row_wise_detection:
    input:
        studies = lambda wc: get_row_wise_input_files(config[wc.dataset], config["id_pattern"], "work_folder/formated/bait_prey_publications.csv"),
    output:
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
    log:
        "logs/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.log",
    threads:
        1
    run:
        studies = list(input.studies)
        study_folder = next((f for f in studies if os.path.isdir(f)), None)
        if study_folder:  # Checkpoint left unresolved (spawned job), select the studies ourselves
            studies = select_row_wise_study_files(
                study_folder, config[wildcards.dataset], config["id_pattern"],
                "work_folder/formated/bait_prey_publications.csv"
            )
        with open(output.row_wise_data, "w") as w:
            w.write("\t".join(["bait", "prey", "experiment", "detection"]) + "\n")
            for study in studies:
                with open(study, "r") as f:
                    next(f) # skip header
                    for l in f:
                        values = l.strip().split("\t")
                        bait = values[0]; prey = values[1]
                        detection = values[3]; exp = f"{values[5]}_{values[4]}"
                        w.write("\t".join([bait, prey, exp, detection]) + "\n")


rule test_distribution:
    params:
        search_space_minimum = [0, 200, 400, 800, 1600],
    input:
        pod_data = "work_folder/analysis/POD/directional/POD_{dataset}.pq",
    output:
        probabilities = "work_folder/analysis/ProteinPromiscuity/detection_distribution/{dataset}_detectability.tsv",
        plot_distributions = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{dataset}_detection_distribution.png",
        plot_bait_prey_cor = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{dataset}_bait_prey_detectability_correlation.png",
    threads:
        20
    script:
        "scripts/plot_detection_distribution.R"


rule test_degree_distribution:
    params:
        pod_limits = [0.02, 0.05, 0.15, 0.3],
        fdr_limits = [0.1, 0.3, 0.5, 0.99],
    input:
        pod_data = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
        fdr_data = "work_folder/analysis/FDR_aware/pair_posterior/{network_type}_{dataset}.pq",
    output:
        pod_degree = "work_folder/analysis/ProteinPromiscuity/detection_distribution/{network_type}_{dataset}_pod_degree.tsv",
        fdr_degree = "work_folder/analysis/ProteinPromiscuity/detection_distribution/{network_type}_{dataset}_fdr_degree.tsv",
        plot_bait_prey_cor_pod = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{network_type}_{dataset}_pod_bait_prey_degree_correlation.png",
        plot_bait_prey_cor_fdr = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{network_type}_{dataset}_fdr_bait_prey_degree_correlation.png",
        plot_distributions_pod = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{network_type}_{dataset}_pod_degree_distribution.png",
        plot_distributions_fdr = "work_folder/analysis/ProteinPromiscuity/detection_distribution/plot/{network_type}_{dataset}_fdr_degree_distribution.png",
    threads:
        20
    script:
        "scripts/plot_degree_distribution.R"