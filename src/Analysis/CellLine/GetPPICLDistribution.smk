import pandas as pd

rule distributions_MS_cl:
    params:
        methods = [
            "MI-0006",
            "MI-0007",
            "MI-1314",
            "MI-0096",
            "MI-0004",
            "MI-0019"
        ],
        min_ppis = 1000
    input:
        ppi_pid_file = "work_folder/intact/method_subset/bait_prey_publications.csv",
        pid_cl_counts = "work_folder/pid_cell_line/pid_ppi_count.csv"
    output:
        single_cl = "work_folder/pid_cl_ppi/single_cl.csv",
        multi_cl = "work_folder/pid_cl_ppi/multi_cl.csv",
        top_multi_pid = "work_folder/pid_cl_ppi/top_multi_cl.csv"
    run:
        ppi_df = pd.read_csv(input.ppi_pid_file, sep="\t")
        ppi_df = ppi_df[ppi_df["detection_method"].isin(params.methods)]
        pids_cl_df = pd.read_csv(input.pid_cl_counts, sep = "\t")

        ppi_pid_cl_df = ppi_df.merge(pids_cl_df, on="pubmed_id")
        single_ss = ppi_pid_cl_df[ppi_pid_cl_df["cl_count"] == 1]
        single_ss.to_csv(output.single_cl, sep="\t", index=False)

        multi_ss = ppi_pid_cl_df[ppi_pid_cl_df["cl_count"] != 1]
        multi_ss.to_csv(output.multi_cl, sep="\t", index=False)

        count_cl = multi_ss.groupby(["pubmed_id", "cl_count"], as_index=False).count()
        count_cl["count"] = count_cl["bait"]
        count_cl_top = count_cl[count_cl["count"] >= 1000]
        count_cl_top[["pubmed_id", "cl_count", "count"]].to_csv(output.top_multi_pid, sep="\t", index=False)





