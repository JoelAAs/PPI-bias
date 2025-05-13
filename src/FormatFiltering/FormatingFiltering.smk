from formating import *

##### Rules
rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab = "data/intact/human.txt"
    output:
        formated = "work_folder/formated/bait_prey_publications.csv"
    run:
        intact_df = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(intact_df)
        bait_prey_df.to_csv(
            output.formated,
            sep="\t",
            index = None
        )

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