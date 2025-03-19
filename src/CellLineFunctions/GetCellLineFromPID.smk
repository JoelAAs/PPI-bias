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
        #pubtator_count["pubmed_id"] = pubtator_count.index
        pubtator_count.reset_index(inplace=True)
        pubtator_count = pubtator_count.rename({"info_type": "cl_count"}, axis = 1)
        pubtator_count.to_csv(
            output.pid_cl_count,
            sep="\t",
            index=False
        )

        pubtator_df = pubtator_df.merge(pubtator_count, on="pubmed_id")
        for cell_line in config["cell_lines"]:
            if cell_line == "all":
                cl_ss_pubtator_df = pubtator_df
            else:
                cl_ss_pubtator_df = pubtator_df[pubtator_df["cl_id"] == cell_line]
            cl_ss_pubtator_df[["pubmed_id", "cl_count"]].to_csv(
                f"work_folder/pid_cell_line/{cell_line}.csv",
                sep="\t",
                index=False
            )



#
# rule get_max_cl_pid_lists:
#     """
#     refactor later can be simplified
#     """
#     input:
#         cl = "work_folder/pid_cell_line/{cell_line}.csv",
#         pid_cl_count = "work_folder/pid_cell_line/pid_ppi_count.csv"
#     output:
#         cl = "work_folder/pid_cell_line//{cell_line}{n}.csv"
#     run:
#         pid_df = pd.read_csv(input.cl, sep = "\t")
#         pid_count_df = pd.read_csv(input.pid_cl_count, sep="\t")
#
#         pid_df = pid_df.merge(pid_count_df, on="pubmed_id", how="inner")
#         pid_df = pid_df[pid_df["cl_count"] <= int(wildcards.n)]
#         print(pid_df)
#         pid_df.to_csv(output.cl, sep="\t", index=False)

