rule coessentiality:
    output:
        gnomad_lof_metrics = "work_folder/analysis/coessentiality/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    shell:
        """
        wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
        -O {output.gnomad_lof_metrics}
        """


rule get_oe_lof_category:
    params:
        low_quantile = config["oe_lof_low_quantile"]
    input:
        gnomad_lof_metrics = "work_folder/analysis/coessentiality/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        symbol_map = "work_folder/analysis/coexpression/uniprot_to_hgnc_symbol.tsv"
    output:
        oe_lof_category = "work_folder/analysis/coessentiality/oe_lof_category.tsv"
    log:
        "logs/analysis/coessentiality/oe_lof_category.log"
    script:
        "scripts/get_oe_lof_category.py"


rule get_balanced_set_coessentiality:
    input:
        max_positive = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        max_negative = f"work_folder/subsets/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv",
        oe_lof_category = "work_folder/analysis/coessentiality/oe_lof_category.tsv"
    output:
        balanced_positive = f"work_folder/analysis/coessentiality/balanced/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        balanced_negative = f"work_folder/analysis/coessentiality/balanced/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv"
    log:
        "logs/analysis/coessentiality/balanced/{dataset}_{network_type}.log"
    script:
        "scripts/balance_coessentiality.py"


rule test_coessentiality_hri_hnri:
    input:
        balanced_positive = f"work_folder/analysis/coessentiality/balanced/{{dataset}}_{{network_type}}_limit_{config['positive_max']}_pos.csv",
        balanced_negative = f"work_folder/analysis/coessentiality/balanced/{{dataset}}_{{network_type}}_limit_{config['negative_max']}_neg.csv",
        oe_lof_category = "work_folder/analysis/coessentiality/oe_lof_category.tsv"
    output:
        coessentiality_OR = "work_folder/analysis/coessentiality/{dataset}_{network_type}_coessentiality_OR.tsv"
    threads: 20
    log:
        "logs/analysis/coessentiality/{dataset}_{network_type}_coessentiality_OR.log"
    script:
        "scripts/test_coessentiality_hri_hnri.py"


rule plot_coessentiality:
    input:
        coessentiality_OR = expand(
            "work_folder/analysis/coessentiality/{dataset}_{{network_type}}_coessentiality_OR.tsv",
            dataset=config["datasets"]
        )
    output:
        "work_folder/analysis/coessentiality/plots/{network_type}_coessentiality_OR.png"
    log:
        "logs/analysis/coessentiality/plots/{network_type}_coessentiality_OR.log"
    script:
        "scripts/plot_coessentiality_or.R"

