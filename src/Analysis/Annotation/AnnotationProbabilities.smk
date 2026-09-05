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


rule protein_annotation_jaccard:
    params:
        id_pattern = config["id_pattern"]
    input:
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        jaccard = "work_folder/analysis/shared_annotation_proportions/jaccard_{dataset}.csv"
    log:
        "logs/analysis/shared_annotation_proportions/jaccard_{dataset}.log"
    script:
        "scripts/protein_annotation_jaccard.py"


rule plot_anotation_jaccard:
    input:
        expand("work_folder/analysis/shared_annotation_proportions/jaccard_{dataset}.csv",
            dataset=config["datasets"])
    output:
        "work_folder/analysis/shared_annotation_proportions/plots/jaccard_mean_vs_pod.png"
    log:
        "logs/analysis/shared_annotation_proportions/plots/jaccard_mean_vs_pod.log"
    script:
        "scripts/plot_annotation_jaccard.py"


rule table_shared_annotations:
    # One row per annotation, with the odds ratio and q-value of the combined (flat), MS and
    # Y2H datasets side by side.
    params:
        dataset_labels = {"flat": "combined", "ms": "ms", "y2h": "y2h"}
    input:
        expand("work_folder/analysis/shared_annotation_proportions/{dataset}_{{network_type}}.tsv",
            dataset=config["datasets"])
    output:
        table = "work_folder/analysis/shared_annotation_proportions/{network_type}_OR_q_table.tsv"
    log:
        "logs/analysis/shared_annotation_proportions/{network_type}_OR_q_table.log"
    run:
        columns = {}
        for path in input:
            dataset = path.split("/")[-1].removesuffix(f"_{wildcards.network_type}.tsv")
            label = params.dataset_labels[dataset]
            df = pd.read_csv(path, sep="\t").set_index(["annotation", "annotation_type"])
            columns[label] = df[["odds_ratio", "q_val"]].rename(
                columns={"odds_ratio": f"OR_{label}", "q_val": f"q_{label}"}
            )

        table = pd.concat(
            [columns[label] for label in params.dataset_labels.values()], axis=1
        ).reset_index()
        table = table.sort_values(["annotation_type", "annotation"])
        # Odds ratios to 3 decimals; q-values in scientific notation, since they bottom out
        # at the permutation floor (~2e-4) and would otherwise all round to 0.000.
        or_columns = [f"OR_{label}" for label in params.dataset_labels.values()]
        q_columns = [f"q_{label}" for label in params.dataset_labels.values()]
        table[or_columns] = table[or_columns].round(3)
        table[q_columns] = table[q_columns].map(lambda v: f"{v:.3e}")
        table.to_csv(output.table, sep="\t", index=False)
