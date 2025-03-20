import pandas as pd

rule get_pid_cla:
    input:
        pubtator = config["pid_context"]
    output:
        cl_subset = expand(
            "work_folder/pid_cell_line/{cell_line}.csv",
            cell_line = config["cell_lines"]
        ),
        pid_cl_count = "work_folder/pid_cell_line/pid_ppi_count.csv"
    run:
        pubtator_df = pd.read_csv(
            input.pubtator,
            sep="\t"
        )
        pubtator_count = pubtator_df.groupby("pubmed_id").count()["info_type"]
        pubtator_count = pd.DataFrame(pubtator_count)
        pubtator_count.reset_index(inplace=True)
        pubtator_count = pubtator_count.rename({"info_type": "cl_count"}, axis = 1)
        pubtator_count.to_csv(
            output.pid_cl_count,
            sep="\t",
            index=False
        )

        pubtator_df = pubtator_df.merge(pubtator_count, on="pubmed_id")
        for cell_line in config["cell_lines"]:
            cl_ss_pubtator_df = pubtator_df[pubtator_df["cl_id"] == cell_line].copy()
            cl_ss_pubtator_df["cell_line"] = cell_line
            cl_ss_pubtator_df[["pubmed_id", "cl_count", "cell_line"]].to_csv(
                f"work_folder/pid_cell_line/{cell_line}.csv",
                sep="\t",
                index=False
            )



