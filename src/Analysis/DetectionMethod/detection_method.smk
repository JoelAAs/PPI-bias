def get_input_files(method, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]

    ppi_df = pd.read_csv(filename, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"] == method]
    ppi_df = ppi_df[
        ppi_df[f"{id_pattern}_bait"] != ppi_df[f"{id_pattern}_prey"]
    ]
    if ppi_df.empty():
        raise ValueError(f"No studies for method: {method}")
    ppi_df = ppi_df[~ppi_df[[f"{id_pattern}_bait", f"{id_pattern}_prey", "pubmed_id"]].duplicated(keep="first")] # Remove isoforms
    if remove_single:
        ppi_df = ppi_df.groupby(["pubmed_id"], as_index=False).size()
        ppi_df = ppi_df[ppi_df["size"] != 1]
    expected = [
        f"{STUDY_FOLDER}/{pid}_{method}.csv" for pid in ppi_df["pubmed_id"].unique()
    ]
    return expected

def multi_method_aggregation(methods):
    return expand(
        rules.aggregate_pids.output.method_aggregate,
        subset=methods
    )

rule aggregate_pids:
    """
    Aggregate data form studies of the same method
    """
    input:
        input_ppi = config["formated_ppi"],
        subsets = lambda wc: multi_method_aggregation(config[wc.subset]) if wc.subset in config else get_input_files(wc.subset, config["id_pattern"], config["formated_ppi"])
    output:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/{subset}_experimental_wise.csv"
    run:
        aggregate_inferred_experiments(input.subsets, output.method_aggregate)
