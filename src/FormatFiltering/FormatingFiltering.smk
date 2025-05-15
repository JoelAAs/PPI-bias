import pandas as pd

from formating import *

##### Rules
rule get_gene_name_uniprot:
    """
    Extract uniprotID-gene_name from miTab aliases column
    """
    input:
        miTab = "data/intact/human.txt"
    output:
        uniprot = "work_folder/intact/uniprot_to_gene_name.csv"
    run:
        get_gene_names(input.miTab, output.uniprot)


rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "data/intact/human.txt",
        gene_names = "work_folder/intact/uniprot_to_gene_name.csv"
    output:
        formated = "work_folder/formated/bait_prey_publications.csv"
    run:
        mitab_df    = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        gene_name_df = pd.read_csv(input.gene_names, sep = "\t")

        mitab_df = mitab_df.merge(gene_name_df, left_on="uniprot_bait", right_on="uniprot_id")
        del mitab_df["uniprot_id"]
        mitab_df = mitab_df.merge(gene_name_df, left_on="uniprot_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del mitab_df["uniprot_id"]

        bait_prey_df.to_csv(
            output.formated,
            sep="\t",
            index = None
        )

