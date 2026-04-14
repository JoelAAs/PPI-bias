import pandas as pd
from formating import *

##### Rules
rule get_intact:
    output:
        intact="work_folder/data/intact/human.txt"
    log:
        "logs/data/intact/human.log"
    shell:
        """
        wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt -O {output.intact} > {log} 2>&1
        """

rule get_gene_name_uniprot:
    """
    Extract uniprotID-gene_name from miTab aliases column
    """
    input:
        miTab = "work_folder/data/intact/human.txt"
    output:
        uniprot = f"work_folder{pn}/gene_names/uniprot_to_gene_name.csv"
    log:
        f"logs{pn}/gene_names/uniprot_to_gene_name.log"
    run:
        get_gene_names(input.miTab, output.uniprot)


rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "work_folder/data/intact/human.txt",
        gene_names_unipriot = f"work_folder{pn}/gene_names/uniprot_to_gene_name.csv",
        gene_names_swissprot= f"work_folder{pn}/gene_names/uniprot_to_sp.csv"
    output:
        formated = f"work_folder{pn}/formated/bait_prey_publications.csv",
        gene_names = f"work_folder{pn}/gene_names/gene_names.csv"
    log:
        f"logs{pn}/formated/bait_prey_publications.log"
    run:
        mitab_df     = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        gene_name_uniprot_df = pd.read_csv(input.gene_names_unipriot, sep = "\t")
        gene_name_swissprot_df = pd.read_csv(input.gene_names_swissprot, sep = "\t")


        gene_name_df = gene_name_uniprot_df.merge(gene_name_swissprot_df,
            left_on="gene_name", right_on="intact_gene_name"
        )
        gene_name_df = gene_name_df[["uniprot_id", "sp_gene_name"]]
        gene_name_df.columns = ["uniprot_id", "gene_name"]
        gene_name_df.to_csv(output.gene_names, sep="\t", index=None)

        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
        del bait_prey_df["uniprot_id"]
        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del bait_prey_df["uniprot_id"]

        bait_prey_df.to_csv(
            output.formated,
            sep="\t",
            index = None
        )