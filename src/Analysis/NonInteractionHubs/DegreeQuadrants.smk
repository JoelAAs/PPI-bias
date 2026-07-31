
rule plot_degree_vs_test:
    input:
        expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        )
    output:
        ref  = "work_folder/analysis/protein_degree/plots/{network_type}_degree_dist.png",
        dist = "work_folder/analysis/protein_degree/plots/{network_type}_degree_vs_pubs.png"
    log:
        "logs/analysis/protein_degree/plots/{network_type}_degree.log"
    script:
        "scripts/plots/plot_degree.R"


rule plot_degree_distributions:
    input:
        degree = expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        ),
        bin_summary = expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree_bin_summary.tsv",
            dataset=config["datasets"]
        )
    output:
        distributions = "work_folder/analysis/protein_degree/plots/{network_type}_degree_distributions.png",
        heatmap = "work_folder/analysis/protein_degree/plots/{network_type}_degree_bin_heatmap.png"
    log:
        "logs/analysis/protein_degree/plots/{network_type}_degree_distributions.log"
    script:
        "scripts/plots/plot_degree_distributions.R"

