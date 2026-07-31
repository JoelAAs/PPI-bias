
rule plot_degree_vs_test:
    input:
        expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        )
    output:
        dist = "work_folder/analysis/protein_degree/plots/{network_type}_degree_vs_pubs.png",
        ref  = "work_folder/analysis/protein_degree/plots/{network_type}_degree_dist.png"
    log:
        "logs/analysis/protein_degree/plots/{network_type}_degree.log"
    script:
        "scripts/plots/plot_degree.R"