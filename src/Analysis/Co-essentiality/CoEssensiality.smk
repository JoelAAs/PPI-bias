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
        max_positive = "work_folder/naive_split/{network_type}/full/pos/{exp_set}.tsv",
        max_negative = "work_folder/naive_split/{network_type}/full/neg/{exp_set}.tsv",
        oe_lof_category = "work_folder/analysis/coessentiality/oe_lof_category.tsv"
    output:
        balanced_positive = "work_folder/analysis/coessentiality/balanced/{exp_set}_{network_type}_pos.csv",
        balanced_negative = "work_folder/analysis/coessentiality/balanced/{exp_set}_{network_type}_neg.csv"
    log:
        "logs/analysis/coessentiality/balanced/{exp_set}_{network_type}.log"
    script:
        "scripts/balance_coessentiality.py"


rule test_coessentiality_hri_hnri:
    input:
        balanced_positive = "work_folder/analysis/coessentiality/balanced/{exp_set}_{network_type}_pos.csv",
        balanced_negative = "work_folder/analysis/coessentiality/balanced/{exp_set}_{network_type}_neg.csv",
        oe_lof_category = "work_folder/analysis/coessentiality/oe_lof_category.tsv"
    output:
        coessentiality_OR = "work_folder/analysis/coessentiality/{exp_set}_{network_type}_coessentiality_OR.tsv"
    threads: 20
    log:
        "logs/analysis/coessentiality/{exp_set}_{network_type}_coessentiality_OR.log"
    script:
        "scripts/test_coessentiality_hri_hnri.py"


rule plot_coessentiality:
    input:
        coessentiality_OR = expand(
            "work_folder/analysis/coessentiality/{exp_set}_{{network_type}}_coessentiality_OR.tsv",
            exp_set=["fdr", "pod"]
        )
    output:
        "work_folder/analysis/coessentiality/plots/{network_type}_coessentiality_OR.png"
    log:
        "logs/analysis/coessentiality/plots/{network_type}_coessentiality_OR.log"
    script:
        "scripts/plot_coessentiality_or.R"

