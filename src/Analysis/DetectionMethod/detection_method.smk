def get_input_files(method, id_pattern, filename, remove_single=True):
    """
    Get list of expected output of single or multiple experiments with negative data inferred
    :param method: (stirng) Either a single MI tag or an aggregation (defined in config)
    :param id_pattern: (string) Either "gene_name" or "uniprot_id"
    :param filename: (string) path to formated bait-prey file
    :param remove_single: (bool) Will remove any bait-prey experiments reporting a single interaction (this does not remove single bait studies)
    :return: (list) list of paths of expected output
    """
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]
    STUDY_FOLDER = "work_folder" + STUDY_FOLDER.split("work_folder")[1]
    ppi_df = pd.read_csv(filename, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"] == method]
    ppi_df = ppi_df[
        ppi_df[f"{id_pattern}_bait"] != ppi_df[f"{id_pattern}_prey"]
    ]
    if ppi_df.empty:
        raise ValueError(f"No studies for method: {method}")
    ppi_df = ppi_df[~ppi_df[[f"{id_pattern}_bait", f"{id_pattern}_prey", "pubmed_id"]].duplicated(keep="first")] # Remove isoforms, if considering gene names
    if remove_single:
        ppi_df = ppi_df.groupby(["pubmed_id"], as_index=False).size()
        ppi_df = ppi_df[ppi_df["size"] != 1]
    expected = [
        f"{STUDY_FOLDER}/{pid}_{method}.csv" for pid in ppi_df["pubmed_id"].unique()
    ]
    return expected

def multi_method_aggregation(methods):
    # hardcoded to play nice with snakemake-fs
    return expand(
        "work_folder/inferred_search_space/aggregated/methods/{subset}_experimental_wise.csv",
        subset=methods
    )

def get_subsets(wc):
    """
    :param wc: snakemake wildcards
    :return: either mulit-method (MS/y2h) aggregation or single method (MI-code) aggregation
    """
    if wc.subset in config:
        return multi_method_aggregation(config[wc.subset])
    return get_input_files(wc.subset,config["id_pattern"],config["formated_ppi"])

rule aggregate_pids:
    """
    Aggregate data from studies of the same method
    """
    params:
        id_pattern = config["id_pattern"]
    input:
        input_ppi = storage.fs(config["formated_ppi"]),
        subsets = lambda wc: get_subsets(wc)
    output:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/{subset}_experimental_wise.csv"
    log:
        "logs/inferred_search_space/aggregated/methods/{subset}_experimental_wise.log"
    run:
        single = wildcards.subset not in config
        aggregate_inferred_experiments(input.subsets, output.method_aggregate, params.id_pattern, single)
