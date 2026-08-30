rule get_other_ppis:
    input:
        miTab = "work_folder/data/intact/human.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        other_ppis = "work_folder/formated/other_method_ppis.csv"
    log:
        "logs/formated/other_method_ppis.log"
    run:
        mitab_df = filter_mitab(input.miTab)
        mitab_df["detection_method"] = mitab_df["detection_method"].str.replace(":", "-")

        selected_methods = set(config["ms"]) | set(config["y2h"])
        other_df = mitab_df[~mitab_df["detection_method"].isin(selected_methods)]

        other_df = other_df[["IDA", "IDB", "detection_method"]].dropna()
        other_df = other_df.rename(columns={"IDA": "prot_a", "IDB": "prot_b"})

        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        other_df = other_df.merge(gene_name_df, left_on="prot_a", right_on="uniprot_id")
        del other_df["uniprot_id"]
        other_df = other_df.merge(gene_name_df, left_on="prot_b", right_on="uniprot_id", suffixes=("_a", "_b"))
        del other_df["uniprot_id"]

        other_df = other_df.drop_duplicates()
        other_df.to_csv(output.other_ppis, sep="\t", index=None)


rule check_FDR:
    input:
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq"
        other_ppi_intact = ""

