rule get_protein_localisation:
    params:
        id_pattern = config["id_pattern"],
    input:
        annotation = config["localisation_file"],
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv",
    output:
        protein_localisation = "work_folder/analysis/gene_names/protein_localisation.tsv",
    log:
        "logs/analysis/gene_names/protein_localisation.log",
    threads:
        1
    script:
        "scripts/get_protein_localisation.py"


rule plot_rank_protein_detectability:
    input:
        bait_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_bait_detectability.tsv",
        prey_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_prey_detectability.tsv",
        protein_localisation = "work_folder/analysis/gene_names/protein_localisation.tsv",
    output:
        plot = "work_folder/analysis/degree_estimation/plot/{dataset}_spearman_cross_detectability_localisation.png",
    threads:
        1
    script:
        "scripts/plot_spearman_cross_detectability.R"


rule plot_bait_prey_detectability_correlation:
    input:
        bait_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_bait_detectability.tsv",
        prey_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_prey_detectability.tsv",
    output:
        plot = "work_folder/analysis/degree_estimation/plot/{dataset}_bait_prey_detectability_correlation.png",
    threads:
        1
    script:
        "scripts/plot_bait_prey_detectability_correlation.R"
