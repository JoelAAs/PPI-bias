
def reformat_long_cl(df, pattern, col_name, id_column="gene_name_prey"):
    df_col = [c for c in df.columns if pattern in c]
    df_long = pd.melt(
        df[[id_column,] + df_col],
        id_vars=id_column,var_name="CVCL",value_name=col_name
    )
    df_long["CVCL"] = df_long["CVCL"].apply(lambda x: x.replace(pattern,""))

    return df_long

def nested_dict():
    return defaultdict(nested_dict)

def get_input_for_aggregation(wc, filename, cell_line_methods):
    ckpt = checkpoints.infer_experimental_search_space.get(cell_line="_cell_line").output[0]
    CL_FOLDER = "work_folder/inferred_search_space/experimental_cell_line" # Explicit since it doesn't work otherwise
    prefix = workflow.storage_settings.default_storage_prefix
    query = f"{prefix.rstrip('/')}/{filename}" if prefix else filename
    storage_file = storage.fs(query)
    storage_object = storage_file.flags["storage_object"]
    storage_object.local_path().parent.mkdir(parents=True, exist_ok=True)  # retrieve_object() rsyncs, it does not create the target dir
    storage_object.retrieve_object()
    cl_df = pd.read_csv(storage_file, sep="\t")
    cl_df = cl_df[
        cl_df[f"gene_name_bait"] != cl_df[f"gene_name_prey"]
        ] # remove bait-bait

    cl_df = cl_df[cl_df["detection_method"].isin(cell_line_methods)]    

    cl_df = cl_df[
        ~cl_df[[
            "gene_name_bait", "gene_name_prey",
            "pubmed_id", "detection_method", "CVCL"
        ]].duplicated(keep="first")] # remove isoforms

    cl_df = cl_df[["pubmed_id", "detection_method", "CVCL"]]
    cl_df = cl_df[cl_df.duplicated(keep=False)]

    expected_input = {
        f"{CL_FOLDER}/{pubmed_id}_{detection_method}_{cl_id}.csv"
        for _, (pubmed_id, detection_method, cl_id) in cl_df.iterrows()
    }
    return expected_input


rule aggregate_inferred_studies_cell_line:
    """
    Aggregate experiments assuming that any prey observed in studies is tested against all baits
    """
    input:
        ppi_file = "work_folder/formated/bait_prey_CVCL.csv",
        cl_pids = lambda wc: get_input_for_aggregation(wc, "work_folder/formated/bait_prey_CVCL.csv", config["cell_line_methods"])
    output:
        cell_line_counts = "work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    log:
        "logs/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.log"
    run:
        aggregate_inferred_experiments(input.cl_pids, output.cell_line_counts, "gene_name", single=True) # must be gene name as bioplex only does genes
