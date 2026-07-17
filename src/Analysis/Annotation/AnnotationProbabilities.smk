import pandas as pd


rule get_annotation_proteins:
    input:
        gene_names  = "work_folder/gene_names/uniprot_to_gene_name.csv",
        annotation  = config["localisation_file"]
    output:
        annotation_proteins = "work_folder/analysis/shared_annotation_proportions/annotation_proteins.csv"
    log:
        "logs/analysis/shared_annotation_proportions/annotation_proteins.log"
    script:
        "scripts/get_annotation_proteins.py"


rule test_shared_annotations:
    params:
        bait_column = f"{config['id_pattern']}_bait",
        prey_column = f"{config['id_pattern']}_prey"
    input:
        gene_names  = "work_folder/gene_names/uniprot_to_gene_name.csv",
        annotation  = "work_folder/analysis/shared_annotation_proportions/annotation_proteins.csv",
        edges_pos = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        edges_neg = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv"
    output:
        annotation_protortions = "work_folder/analysis/shared_annotation_proportions/{dataset}_{network_type}.tsv"
    threads: 20
    log:
        "logs/analysis/shared_annotation_proportions/{dataset}_{network_type}.log"
    script:
        "scripts/test_shared_annotations.py"


rule plot_shared_annotations:
    input:
        expand("work_folder/analysis/shared_annotation_proportions/{dataset}_{{network_type}}.tsv",
            dataset=config["datasets"])
    output:
        "work_folder/analysis/shared_annotation_proportions/plots/{network_type}_OR.png"
    log:
        "logs/analysis/shared_annotation_proportions/plots/{network_type}_OR.log"
    script:
        "scripts/plot_shared_annotations.R"
