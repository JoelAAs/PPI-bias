from localisation_support import *
import numpy as np


rule compare_HCI_negatome:
    params:
        localisation_csv = config["localisation_file"]
    input:
        experimental_negatome = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv",
        hci = "work_folder/inferred_search_space/analysis/bias_reduced_ppis/high_confidence.csv"
    output:
        test_results = "work_folder/inferred_search_space/analysis/localisation/HCL_vs_negatome_table.csv",
        test_table = "work_folder/inferred_search_space/analysis/localisation/HCL_vs_negatome_test.txt"
    run:
        df_negatome  = pd.read_csv(input.experimental_negatome,  sep="\t")
        df_negatome["group"] = "negatome"
        df_hci = pd.read_csv(input.hci, sep="\t")
        df_hci["group"] = "hci"

        df_localisation = pd.read_csv(params.localisation_csv, sep="\t")

        df_negatome = add_localisation(df_negatome, df_localisation)
        df_hci = add_localisation(df_hci, df_localisation)
        full_df = pd.concat([df_negatome, df_hci])

        test_df = full_df.groupby(
            ["group", "localisation_match"], as_index=False
        ).size().sort_values(by="localisation_match", ascending= False).sort_values(by="group")
        test_table = np.reshape(test_df["size"], (2,2))
        loc_or, loc_pval = fisher_exact(test_table)

        test_df.to_csv(output.test_table, sep="\t", index=False)
        with open(output.test_results, "w") as w:
            w.write("odds_ratio\tpvalue\n")
            w.write(f"{loc_or}\t{loc_pval}\n")





