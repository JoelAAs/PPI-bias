
rule get_pid_cla:
    input:
        pubtator = config["pid_context"]
    output:
        expand(
            "work_folder/pid_cell_line/{cell_line}.csv",
            cell_line = config["cell_lines"]
        )
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



