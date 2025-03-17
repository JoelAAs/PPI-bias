import re
import numpy as np
import pandas as pd


def get_overlap_input(wc):
    expected_input = [
        f"work_folder/intact/pair_count/ppi_pair_counts_{c_cl}_{method}.csv"
        for method in config[wc.detection_method]
        for c_cl in config["cell_lines"]
        if c_cl != "all"
    ]
    return expected_input



rule get_cell_line_ppi_overlap:
    params:
        bait_prey_specific = True
    input:
        ppi_subset = lambda wc: get_overlap_input(wc)
    output:
        iou = "work_folder/ppi_cl_overlap/IoU_{detection_method}.csv",
        intersect_n = "work_folder/ppi_cl_overlap/intersection_{detection_method}.csv",
        unique_n = "work_folder/ppi_cl_overlap/unique_{detection_method}.csv",
        one_hot = "work_folder/ppi_cl_overlap/one_hot_{detection_method}.csv"

    run:
        cl_ppi_directional_dict = dict()
        all_directional = set()
        cl_ppi_dict = dict()
        all_undirectional = set()

        for cl_ppi in input.ppi_subset:
            cl = re.search(
                r"ppi_pair_counts_([a-zA-Z0-9-_]+)_[A-Z0-9-]+.csv",
                cl_ppi).groups()[0]
            cl_df = pd.read_csv(
                cl_ppi,
                sep="\t"
            )
            if not cl in  cl_ppi_directional_dict:
                cl_ppi_directional_dict[cl] = set()
                cl_ppi_dict[cl] = set()

            cl_df = cl_df[cl_df["observed_interactions"] != 0]
            if not cl_df.empty:
                cl_ppi_directional_dict[cl].update(set(
                    cl_df[["bait", "prey"]].apply("|".join, axis=1).tolist()
                ))
                all_directional.update(cl_ppi_directional_dict[cl])

                sorted_id = lambda x: "|".join(sorted(x))
                cl_ppi_dict[cl].update(set(
                    cl_df[["bait", "prey"]].apply(sorted_id, axis=1).tolist()
                ))
                all_undirectional.update(cl_ppi_dict[cl])

        cls = list(cl_ppi_dict.keys())
        iou_matrix = np.zeros((len(cls), len(cls)))
        n_shared_matrix = np.zeros((len(cls), len(cls)))
        n_unique_matrix = np.zeros((len(cls), len(cls)))

        for i in range(len(cls)-1):
            for j in range(i + 1, len(cls)):
                intersect_ppi = cl_ppi_dict[cls[i]].intersection(
                    cl_ppi_dict[cls[j]])
                union_ppi = cl_ppi_dict[cls[i]].union(cl_ppi_dict[cls[j]])

                iou_matrix[i, j] = len(intersect_ppi)/len(union_ppi)
                iou_matrix[j, i] = iou_matrix[i, j]

                n_shared_matrix[i, j] = len(intersect_ppi)
                n_shared_matrix[j, i] = n_shared_matrix [i, j]

                n_unique_matrix[i, j] = len(cl_ppi_dict[cls[i]] - intersect_ppi)
                n_unique_matrix[j, i] = len(cl_ppi_dict[cls[j]] - intersect_ppi)

        #diagonal
        for i in range(len(cls)):
            iou_matrix[i, i] = 1
            n_shared_matrix[i, i] = len(cl_ppi_dict[cls[i]])
            n_unique_matrix[i, i] = 0


        pd.DataFrame(
            iou_matrix, index=cls, columns=cls).to_csv(
            output.iou, sep="\t")
        pd.DataFrame(
            n_shared_matrix, index=cls, columns=cls).to_csv(
            output.intersect_n, sep="\t")
        pd.DataFrame(
            n_unique_matrix, index=cls, columns=cls).to_csv(
            output.unique_n, sep="\t")

        one_hot_rows = []
        for cl in cls:
            row = {ppi: 1 for ppi in cl_ppi_dict[cl]}
            row.update({"cell_line": cl})
            one_hot_rows.append(row)
        one_hot_df = pd.DataFrame(one_hot_rows)
        one_hot_df.fillna(0, inplace=True)
        one_hot_df.to_csv(output.one_hot, sep="\t", index=False)