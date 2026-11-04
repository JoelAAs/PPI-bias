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
        keep_non_canonical=config["keep_non_canonical"]
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
            params.keep_non_canonical
        )


rule format_miTab:
    """
    Filter and format miTab interaction file into bait-prey-publication-detection_method csv
    """
    params:
        pseudo_replicated = config["pid_pseudo_replication"],
        selected_detection_methods = config["ms"] + config["y2h"]
    input:
        miTab      = "work_folder/data/intact/human.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        formated   = "work_folder/formated/bait_prey_publications.csv"
    log:
        "logs/formated/bait_prey_publications.log"
    run:
        mitab_df     = filter_mitab(input.miTab)
        bait_prey_df = reform_to_bait_prey(mitab_df)
        bait_prey_df = bait_prey_df[~bait_prey_df.pubmed_id.isin(params.pseudo_replicated)]
        bait_prey_df = flag_and_drop_duplicate_studies(bait_prey_df, params.selected_detection_methods)

        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
        del bait_prey_df["uniprot_id"]
        bait_prey_df = bait_prey_df.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del bait_prey_df["uniprot_id"]

        bait_prey_df.to_csv(output.formated, sep="\t", index=None)


                        
def flag_and_drop_duplicate_studies(bp_df, selected_detection_methods, max_jaccard=0.95):
    """
    Some studies only report the intersection between two methods, keep only one.
    Kept in order:
        methods in selected_detection_methods
        then most interactions
        then by removal_order (earliest removed first).
    """
    removal_order = ["MI-1356", "MI-0397", "MI-2215", "MI-0398", "MI-1112", "MI-1113"] # prioriize keeping MS methods, no overlap among MS anymore
    rank = lambda m: removal_order.index(m) if m in removal_order else len(removal_order)

    to_remove = []
    for pid, pid_df in bp_df.groupby("pubmed_id"):
        pairs = {
            m: set(zip(m_df.uniprot_id_bait, m_df.uniprot_id_prey))
            for m, m_df in pid_df.groupby("detection_method")
        }
        kept = []
        for m in sorted(pairs, key=lambda m: (m not in selected_detection_methods, -len(pairs[m]), -rank(m), m)):
            if any(len(pairs[m] & pairs[k]) / len(pairs[m] | pairs[k]) > max_jaccard for k in kept):
                to_remove.append((pid, m))
            else:
                kept.append(m)

    pid_method = pd.MultiIndex.from_frame(bp_df[["pubmed_id", "detection_method"]])
    return bp_df[~pid_method.isin(to_remove)]