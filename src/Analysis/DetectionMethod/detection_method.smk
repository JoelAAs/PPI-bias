def get_input_files(method, filename, remove_single=True):
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get().output[0]

    ppi_df = pd.read_csv(filename, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"] == method]
    if remove_single:
        ppi_df = ppi_df[ppi_df["pubmed_id"].duplicated()] # NOTE: Should this add to method too?
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
        method = lambda wc: get_input_files(wc.single_method, config["ppi_df"])
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