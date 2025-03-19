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
        pub_dict = dict()
        with open(input.pubtator, "r") as f:
            for line in f:
                values = line.strip().split("\t")
                cell_line = values[2]
                pid = values[0]
                if cell_line in pub_dict:
                    pub_dict[cell_line].append(pid)
                else:
                    pub_dict[cell_line] = [pid, ]

        for cell_line in config["cell_lines"]:
            print(cell_line)
            with open(f"work_folder/pid_cell_line/{cell_line}.csv", "w") as w:
                w.write("pubmed_id\n")
                if cell_line == "all":
                    continue

                for pid in pub_dict[cell_line]:
                    w.write(f"{pid}\n")

        pubtator_df = pd.read_csv(
            input.pubtator,
            sep="\t"
        )
        pubtator_count = pubtator_df.groupby("pubmed_id").count()["info_type"]
        pubtator_count = pd.DataFrame(pubtator_count)
        pubtator_count["pubmed_id"] = pubtator_count.index
        pubtator_count = pubtator_count.rename({"info_type": "cl_count"}, axis = 1)
        pubtator_count.to_csv(
            output.pid_cl_count,
            sep="\t",
            index=False
        )


