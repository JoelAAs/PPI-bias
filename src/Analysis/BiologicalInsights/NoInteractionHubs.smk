rule define_hub_set:
    input:
        degree="work_folder/analysis/protein_degree/{dataset}_{network_type}_degree.pq"
    output:
        hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv",
        reference_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_reference_set.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.log"
    script:
        "scripts/define_hub_set.py"


rule fetch_disorder_content:
    input:
        hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv",
        reference_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_reference_set.tsv",
        uniprot_annotation=config["uniprot_annotation"]
    output:
        disorder="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.log"
    script:
        "scripts/fetch_disorder_content.py"


rule test_hub_properties:
    params:
        bootstrap_replicates=config["bootstrap_replicates"],
        seed=config["seed"]
    input:
        hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv",
        reference_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_reference_set.tsv",
        uniprot_annotation=config["uniprot_annotation"],
        disorder="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.tsv"
    output:
        enrichment="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_enrichment.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_enrichment.log"
    script:
        "scripts/test_hub_properties.py"


rule plot_hub_properties:
    input:
        expand(
            "work_folder/analysis/no_interaction_hubs/{dataset}_{{network_type}}_enrichment.tsv",
            dataset=config["datasets"]
        )
    output:
        "work_folder/analysis/no_interaction_hubs/plots/{network_type}_hub_properties.png"
    log:
        "logs/analysis/no_interaction_hubs/plots/{network_type}_hub_properties.log"
    script:
        "scripts/plot_hub_properties.R"
