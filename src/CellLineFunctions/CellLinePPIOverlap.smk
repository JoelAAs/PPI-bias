import pandas as pd

def get_methods_input(wc):
    method_files = [
        f"work_folder/intact/method_subset/{method}.csv"
        for method in config[wc.detection_method]]
    return method_files



rule get_observed_cell_line_ppi_overlap:
    params:
        bait_prey_specific = True
    input:
        pid_cl_list = expand(
            "work_folder/pid_cell_line/{cell_line}.csv", cell_line = config["cell_lines"]
        ),
        method_subsets = lambda wc: get_methods_input(wc)

    output:
        shared_ppi = "work_folder/method_cl_overlap/shared_ppi_{detection_method}.csv",
        unique_ppi = "work_folder/method_cl_overlap/unique_ppi_{detection_method}.csv",
        iou_ppi = "work_folder/method_cl_overlap/iou_ppi_{detection_method}.csv",
        total_ppi = "work_folder/method_cl_overlap/total_ppi_{detection_method}.csv",
        shared_pid = "work_folder/method_cl_overlap/shared_pid_{detection_method}.csv",
        unique_pid = "work_folder/method_cl_overlap/unique_pid_{detection_method}.csv",
        iou_pid= "work_folder/method_cl_overlap/iou_pid_{detection_method}.csv"

    run:
        df_list = []
        for method_ppi in input.method_subsets:
            for pid_cl in input.pid_cl_list:
                pid_df = pd.read_csv(pid_cl,
                    sep="\t",
                    dtype={"pubmed_id": str, "cl_count": int}
                )
                method_df = pd.read_csv(
                    method_ppi,
                    sep="\t",
                    dtype="str"
                )
                single_method_cl_df = method_df.merge(pid_df,on="pubmed_id",how="inner")
                df_list.append(single_method_cl_df)

        methods_cl_df = pd.concat(df_list)

        sorted_id = lambda x: "|".join(sorted(x))
        methods_cl_df["ppi_id"] = methods_cl_df[["bait", "prey"]].apply(sorted_id,axis=1)
        methods_cl_df.to_csv("test.csv", sep ="\t")
        print(methods_cl_df["cl_count"])
        print("här")
        outputs_filenames = {
            "shared_ppi": output.shared_ppi,
            "unique_ppi": output.unique_ppi,
            "iou_ppi": output.iou_ppi,
            "total_ppi": output.total_ppi,
            "shared_pid": output.shared_pid,
            "unique_pid": output.unique_pid,
            "iou_pid": output.iou_pid
        }
        file_dict = {
            key: open(value, "w") for key, value in outputs_filenames.items()
        }

        for max_pid_cl in methods_cl_df["cl_count"].unique():
            cl_count_max_ss = methods_cl_df[methods_cl_df["cl_count"] <= max_pid_cl]
            # ppi overlap
            uniq_ppi_cl = cl_count_max_ss.groupby(["cell_line"])["ppi_id"].unique()
            count_ppi_cl = cl_count_max_ss.groupby(["cell_line"])["ppi_id"].count()

            uniq_pid_cl = cl_count_max_ss.groupby(["cell_line"])["pubmed_id"].unique()

            cls = uniq_ppi_cl.index

            # some symmetric some not
            for cl_from in cls:
                file_dict["total_ppi"].write(f"{cl_from}\t{count_ppi_cl[cl_from]}\t{max_pid_cl}\n")
                for cl_to in cls:
                    uniq_cl_from = set(uniq_ppi_cl[cl_from])
                    uniq_cl_to = set(uniq_ppi_cl[cl_to])
                    intersect_ppi_cls = uniq_cl_from.intersection(uniq_cl_to)
                    union_ppi_cls = uniq_cl_from.union(uniq_cl_to)
                    iou_ppi_cls = len(intersect_ppi_cls)/len(union_ppi_cls)

                    file_dict["shared_ppi"].write(f"{cl_from}\t{cl_to}\t{len(intersect_ppi_cls)}\t{max_pid_cl}\n")
                    file_dict["unique_ppi"].write(
                        f"{cl_from}\t{cl_to}\t{len(uniq_cl_from - intersect_ppi_cls)}\t{max_pid_cl}\n")
                    file_dict["iou_ppi"].write(f"{cl_from}\t{cl_to}\t{iou_ppi_cls}\t{max_pid_cl}\n")

                    uniq_pid_from = set(uniq_pid_cl[cl_from])
                    uniq_pid_to = set(uniq_pid_cl[cl_to])
                    intersect_pid_cls = uniq_pid_from.intersection(uniq_pid_to)
                    union_pid_cls = uniq_pid_from.union(uniq_pid_to)
                    iou_pid_cls = len(intersect_pid_cls)/len(union_pid_cls)
                    file_dict["shared_pid"].write(f"{cl_from}\t{cl_to}\t{len(intersect_pid_cls)}\t{max_pid_cl}\n")
                    file_dict["unique_pid"].write(
                        f"{cl_from}\t{cl_to}\t{len(uniq_pid_from - intersect_pid_cls)}\t{max_pid_cl}\n")
                    file_dict["iou_pid"].write(f"{cl_from}\t{cl_to}\t{iou_pid_cls}\t{max_pid_cl}\n")

        for key, value in file_dict.items():
            value.close()
