ARCHS4_CORRELATION_PKL_URL = "https://s3.amazonaws.com/mssm-data/human_correlation_v2.4.pkl"


rule download_archs4_correlation:
    output:
        pkl="work_folder/data/archs4/human_correlation_v2.4.pkl"
    log:
        "logs/data/archs4/download_correlation.log"
    retries: 3
    shell:
        """
        mkdir -p $(dirname {output.pkl})
        wget "{ARCHS4_CORRELATION_PKL_URL}" -O {output.pkl} > {log} 2>&1
        """


rule map_uniprot_to_hgnc_symbol:
    """
    ARCHS4 keys on HGNC symbol; work_folder/gene_names/uniprot_to_gene_name.csv's
    "gene_name" column is Entrez ID (see script docstring)
    """
    input:
        gene_names="work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        symbol_map="work_folder/analysis/coexpression/uniprot_to_hgnc_symbol.tsv"
    log:
        "logs/analysis/coexpression/uniprot_to_hgnc_symbol.log"
    script:
        "scripts/map_uniprot_to_hgnc_symbol.py"


rule build_coexpression:
    input:
        pkl="work_folder/data/archs4/human_correlation_v2.4.pkl",
        symbol_map="work_folder/analysis/coexpression/uniprot_to_hgnc_symbol.tsv"
    output:
        corr_npy="work_folder/analysis/coexpression/coexpr_matrix.npy",
        gene_index="work_folder/analysis/coexpression/coexpr_gene_index.tsv",
        metadata="work_folder/analysis/coexpression/coexpr_metadata.tsv"
    resources:
        mem_gb=32
    log:
        "logs/analysis/coexpression/build_coexpression.log"
    script:
        "scripts/build_coexpression.py"


rule get_balanced_set:
    input:
        max_positive = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        max_negative = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv",
        gene_index="work_folder/analysis/coexpression/coexpr_gene_index.tsv",
    output:
        balanced_positive = f"work_folder/expression/balanced/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        balanced_negative = f"work_folder/expression/balanced/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv"
    script:
        "scripts/balance_coexpression.py"
        

rule test_hri_vs_hrni:
    input:
        balanced_positive = f"work_folder/expression/balanced/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        balanced_negative = f"work_folder/expression/balanced/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv",
        corr_npy="work_folder/analysis/coexpression/coexpr_matrix.npy",
        gene_index="work_folder/analysis/coexpression/coexpr_gene_index.tsv"
    output:
        hrni_vs_hri="work_folder/analysis/coexpression/{dataset}_{network_type}_hri_hnri.tsv"
    threads: 10
    log:
        "logs/analysis/coexpression/{dataset}_{network_type}.log"
    script:
        "scripts/test_coexpression_hri_hnri.py"




rule protein_mean_coexpression:
    """
    Per-protein mean co-expression with its tested partners, over the full POD table
    (same edge population as protein_degrees.py's sum_tests) - see script docstring.
    """
    params:
        bait_col=f"{config['id_pattern']}_bait",
        prey_col=f"{config['id_pattern']}_prey"
    input:
        pod="work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
        corr_npy="work_folder/analysis/coexpression/coexpr_matrix.npy",
        gene_index="work_folder/analysis/coexpression/coexpr_gene_index.tsv"
    output:
        protein_coexpr="work_folder/analysis/coexpression/{dataset}_{network_type}_protein_mean_coexpr.tsv"
    log:
        "logs/analysis/coexpression/{dataset}_{network_type}_protein_mean_coexpr.log"
    script:
        "scripts/protein_coexpression.py"


rule plot_coexpression:
    input:
        pairs=expand(
            "work_folder/analysis/coexpression/{dataset}_{{network_type}}_pairs.tsv",
            dataset=config["datasets"]
        ),
        degree=expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        ),
        protein_coexpr=expand(
            "work_folder/analysis/coexpression/{dataset}_{{network_type}}_protein_mean_coexpr.tsv",
            dataset=config["datasets"]
        )
    output:
        density="work_folder/analysis/coexpression/plots/{network_type}_coexpression.png",
        sum_tests_vs_coexpression="work_folder/analysis/coexpression/plots/{network_type}_sum_tests_vs_coexpression.png"
    log:
        "logs/analysis/coexpression/plots/{network_type}_coexpression.log"
    script:
        "scripts/plot_coexpression.R"


rule plot_coexpression_stats:
    """
    Summary-level plots (one point per dataset, not per-pair densities):
      - odds ratio of "high" co-expression membership, HRI vs HRNI
      - HRNI summed co-expression vs its uniform-random-pair null
    """
    input:
        summary=expand(
            "work_folder/analysis/coexpression/{dataset}_{{network_type}}_summary.tsv",
            dataset=config["datasets"]
        )
    output:
        or_high="work_folder/analysis/coexpression/plots/{network_type}_or_high_membership.png",
        summed="work_folder/analysis/coexpression/plots/{network_type}_summed_coexpression_vs_random.png"
    log:
        "logs/analysis/coexpression/plots/{network_type}_coexpression_stats.log"
    script:
        "scripts/plot_coexpression_stats.R"
