
rule test_attention_confound:
    params:
        seed=config["seed"],
        bootstrap_replicates=config["bootstrap_replicates"]
    input:
        degree="work_folder/analysis/protein_degree/{dataset}_{network_type}_degree.pq"
    output:
        attention="work_folder/analysis/degree_quadrants/{dataset}_{network_type}_attention.tsv"
    log:
        "logs/analysis/degree_quadrants/{dataset}_{network_type}_attention.log"
    script:
        "scripts/test_attention_confound.py"


rule go_enrichment_quadrants:
    input:
        "work_folder/analysis/degree_quadrants/{dataset}_{network_type}_bins.tsv"
    output:
        "work_folder/analysis/degree_quadrants/{dataset}_{network_type}_go_enrichment.tsv"
    log:
        "logs/analysis/degree_quadrants/{dataset}_{network_type}_go_enrichment.log"
    conda:
        "do_enrichment"
    script:
        "scripts/go_enrichment_quadrants.R"


rule plot_quadrants:
    input:
        expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        )
    output:
        "work_folder/analysis/degree_quadrants/plots/{network_type}_quadrants.png"
    log:
        "logs/analysis/degree_quadrants/plots/{network_type}_quadrants.log"
    script:
        "scripts/plots/plot_quadrants.R"


rule plot_degree_diff_vs_references:
    input:
        expand(
            "work_folder/analysis/degree_quadrants/{dataset}_{{network_type}}_bins.tsv",
            dataset=config["datasets"]
        )
    output:
        diff_vs_references="work_folder/analysis/degree_quadrants/plots/{network_type}_degree_diff_vs_references.png",
        degree_vs_tests="work_folder/analysis/degree_quadrants/plots/{network_type}_degree_vs_tests.png",
        sum_tests_density="work_folder/analysis/degree_quadrants/plots/{network_type}_sum_tests_density.png",
        test_vs_degree="work_folder/analysis/degree_quadrants/plots/{network_type}_sum_tests_vs_degree.png"
    log:
        "logs/analysis/degree_quadrants/plots/{network_type}_degree_diff_vs_references.log"
    script:
        "scripts/plots/plot_degree_diff_vs_references.R"



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
