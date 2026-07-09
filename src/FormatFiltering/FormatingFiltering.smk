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

rule get_sec_ac_mapping:
    output: "work_folder/data/uniprot/sec_ac.txt"
    log: "logs/data/uniprot/sec_ac.log"
    shell:
        """
        wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt -O {output} > {log} 2>&1
        """

rule get_uniprot_idmapping:
    output: "work_folder/data/uniprot/HUMAN_9606_idmapping.dat.gz"
    log: "logs/data/uniprot/HUMAN_9606_idmapping.log"
    shell:
        """
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz -O {output} > {log} 2>&1
        """

rule get_ncbi_gene_info:
    output: "work_folder/data/ncbi/Homo_sapiens.gene_info.gz"
    log: "logs/data/ncbi/Homo_sapiens.gene_info.log"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz -O {output} > {log} 2>&1
        """

rule map_uniprot_to_entrez:
    params:
        drop_isoforms = config["drop_isoforms"]
    input:
        miTab     = "work_folder/data/intact/human.txt",
        sec_ac    = "work_folder/data/uniprot/sec_ac.txt",
        idmapping = "work_folder/data/uniprot/HUMAN_9606_idmapping.dat.gz",
        gene_info = "work_folder/data/ncbi/Homo_sapiens.gene_info.gz"
    output:
        uniprot   = "work_folder/gene_names/uniprot_to_gene_name.csv",
        unmapped  = "work_folder/gene_names/uniprot_unmapped.csv"
    log:
        "logs/gene_names/uniprot_to_gene_name.log"
    run:
        build_entrez_mapping(
            input.miTab, input.sec_ac,
            input.idmapping, input.gene_info,
            output.uniprot, output.unmapped,
            keep_non_canonical = not params.drop_isoforms
        )


rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    input:
        miTab      = "work_folder/data/intact/human.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        formated   = "work_folder/formated/bait_prey_publications.csv",
        gene_names = "work_folder/gene_names/gene_names.csv"
    log:
        "logs/formated/bait_prey_publications.log"
    run:
        mitab_df     = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        gene_name_df.to_csv(output.gene_names, sep="\t", index=None)

        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
        del bait_prey_df["uniprot_id"]
        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del bait_prey_df["uniprot_id"]

        bait_prey_df.to_csv(output.formated, sep="\t", index=None)
