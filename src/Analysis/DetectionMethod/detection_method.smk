def get_input_files(method, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]

    ppi_df = pd.read_csv(filename, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"] == method]
    ppi_df = ppi_df[
        ppi_df[f"{id_pattern}_bait"] != ppi_df[f"{id_pattern}_prey"]
    ]
    ppi_df = ppi_df[~[f"{id_pattern}_bait", f"{id_pattern}_prey", "pubmed_id"].duplicated(keep="first")] # Remove isoforms
    if remove_single:
        ppi_df = ppi_df.groupby(["pubmed_id"], as_index=False).size()
        ppi_df = ppi_df[ppi_df["size"] != 1]
    expected = [
        f"{STUDY_FOLDER}/{pid}_{method}.csv" for pid in ppi_df["pubmed_id"].unique()
    ]
    return expected

def multi_method_aggregation(methods):
    return expand(
        rules.aggregate_single_method.output.method_aggregate,
        single_method=methods
    )

rule aggregate_single_method:
    """
    Aggregate data form studies of the same method
    """
    input:
        input_ppi = config["formated_ppi"],
        method = lambda wc: get_input_files(wc.single_method, config["id_pattern"], config["formated_ppi"])
    output:
        method_aggregate = "work_folder/inferred_search_space/aggregated/methods/{single_method}_experimental_wise.csv"
    run:
        aggregate_inferred_experiments(input.method, output.method_aggregate)

rule aggregate_methods:
    """
    Aggregate from groups of methods
    """
    input:
        single_aggregate = lambda wc: multi_method_aggregation(config[wc.multi_method])
    output:
        multi_method = "work_folder/inferred_search_space/aggregated/multi_methods/{multi_method}_experimental_wise.csv"
    run:
        aggregate_inferred_experiments(input.single_aggregate,output.multi_method)