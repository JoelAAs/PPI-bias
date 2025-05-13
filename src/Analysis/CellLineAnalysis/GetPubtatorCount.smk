import re
import pandas as pd


## refactor

rule get_min_prot_cl_count:
    input:
        pid_cl_count = "work_folder/pid_cell_line/pid_ppi_count.csv",
        cl_ppi_pid = "work_folder/intact/method_subset/cell_line_subset/interactions_{cell_line}_{method}.csv"
    output:
        baits = "work_folder/intact/pubtator_cl_count/interaction_min_{cell_line}_{method}_baits.csv",
        preys = "work_folder/intact/pubtator_cl_count/interaction_min_{cell_line}_{method}_preys.csv"
    run:
        ppi_df = pd.read_csv(input.cl_ppi_pid, sep = "\t")
        pid_cl_count_df = pd.read_csv(input.pid_cl_count, sep = "\t")

        ppi_df = ppi_df.merge(pid_cl_count_df, on = "pubmed_id")
        min_bait_dict = dict()
        min_prey_dict = dict()

        for i, row in ppi_df.iterrows():
            match_list = [
                [row["bait"], row["prey"]],
                [row["cl_count"], row["cl_count"]],
                [min_bait_dict, min_prey_dict]
            ]
            for protein, count, prot_dict in zip(*match_list):
                try:
                    if count < prot_dict[protein]:
                        prot_dict[protein] = count
                except KeyError:
                    prot_dict[protein] = count

        for min_prot_dict, outputfile in zip([min_bait_dict, min_prey_dict], output):
            with open(outputfile, "w") as w:
                w.write("uniprot_id\tmin_cl_count\n")
                for protein, value in min_prot_dict.items():
                    w.write(f"{protein}\t{value}\n")

rule percent_unique:
    input:
        min_cl_interactions = expand(
            "work_folder/intact/pubtator_cl_count/interaction_min_{cell_line}_{{method}}_{type}.csv",
            cell_line = config["cell_lines"], type=["baits", "preys"]),
    output:
        overlap = "work_folder/intact/pubtator_cl_count/overlap/{method}.csv"
    run:
        full_cl_count = pd.DataFrame()
        for input_file in input.min_cl_interactions:
            cl_pattern = f"work_folder/intact/pubtator_cl_count/interaction_min_([a-zA-Z0-9-_]+)_{wildcards.method}_"
            cl = re.search(cl_pattern, input_file).groups()[0]
            type_pattern = f"work_folder/intact/pubtator_cl_count/interaction_min_{cl}_{wildcards.method}_([a-z]+).csv"
            type = re.search(type_pattern, input_file).groups()[0]
            if cl != "all":
                cl_count_df = pd.read_csv(input_file, sep="\t")
                cl_count_df["cl"] = cl
                cl_count_df["type"] = type
                full_cl_count = pd.concat([full_cl_count, cl_count_df])

        with open(output.overlap, "w") as w:
            w.write(
                "\t".join([
                    "cell_line",
                    "max_cl_count",
                    "type",
                    "n_intersection",
                    "n_unique",
                ]) + "\n"
            )

            full_cl_count.to_csv("test.csv", sep = "\t", index=False)

            for type in full_cl_count["type"].unique():
                ss_type_df = full_cl_count[full_cl_count["type"] == type]

                for cl_count in ss_type_df["min_cl_count"].unique():
                    ss_count_df = ss_type_df[ss_type_df["min_cl_count"] <= cl_count]

                    for cl in ss_count_df["cl"].unique():
                        ss_count_rest_df = ss_count_df[ss_count_df["cl"] != cl]
                        rest_proteins = set(ss_count_rest_df["uniprot_id"].tolist())
                        ss_count_cl_df = ss_count_df[ss_count_df["cl"] == cl]
                        cl_proteins = set(ss_count_cl_df["uniprot_id"].tolist())


                        n_intersect_protein = len(rest_proteins.intersection(cl_proteins))

                        n_unique = len(cl_proteins) - n_intersect_protein

                        w.write(
                            "\t".join([
                                cl,
                                str(cl_count),
                                type,
                                str(n_intersect_protein),
                                str(n_unique)
                            ]) + "\n"
                        )




