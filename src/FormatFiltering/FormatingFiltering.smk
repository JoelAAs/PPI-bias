from formating import *

##### Rules
rule get_intact:
    output:
        intact="work_folder/data/intact/human.txt"
    shell:
        """
        wget https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt -O {output.intact}
        """

rule get_gene_name_uniprot:
    """
    Extract uniprotID-gene_name from miTab aliases column
    """
    input:
        miTab = "work_folder/data/intact/human.txt"
    output:
        uniprot = f"work_folder{pn}/intact/uniprot_to_gene_name.csv"
    run:
        get_gene_names(input.miTab, output.uniprot)


rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "work_folder/data/intact/human.txt",
        gene_names = f"work_folder{pn}/intact/uniprot_to_gene_name.csv"
    output:
        formated = f"work_folder{pn}/formated/bait_prey_publications.csv"
    run:
        mitab_df    = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        gene_name_df = pd.read_csv(input.gene_names, sep = "\t")

        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
        del bait_prey_df["uniprot_id"]
        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del bait_prey_df["uniprot_id"]

        bait_prey_df.to_csv(
            output.formated,
            sep="\t",
            index = None
        )

