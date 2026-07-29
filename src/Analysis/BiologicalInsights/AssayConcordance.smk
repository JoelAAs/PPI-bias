# Module 4 - is agreement between Y2H and MS concentrated on negatives? Weak/transient
# and co-complex interactions are both real, so positives are where the two assays
# differ in sensitivity, whereas a genuine non-interaction has nothing for either
# assay to detect. See scripts/test_concordance.py docstring for methodology.

wildcard_constraints:
    positive_definition="(any_pos|strict)"


rule build_joint_space:
    params:
        bait_col=f"{config['id_pattern']}_bait",
        prey_col=f"{config['id_pattern']}_prey",
        min_tested=config["concordance_min_tested"]
    input:
        pod_ms="work_folder/analysis/POD/undirectional/POD_ms.pq",
        pod_y2h="work_folder/analysis/POD/undirectional/POD_y2h.pq"
    output:
        joint_space="work_folder/analysis/assay_concordance/joint_space.pq"
    log:
        "logs/analysis/assay_concordance/joint_space.log"
    script:
        "scripts/build_joint_space.py"


rule test_concordance:
    params:
        config=config
    input:
        joint_space="work_folder/analysis/assay_concordance/joint_space.pq"
    output:
        concordance="work_folder/analysis/assay_concordance/concordance_{positive_definition}.tsv"
    threads: 15
    log:
        "logs/analysis/assay_concordance/concordance_{positive_definition}.log"
    script:
        "scripts/test_concordance.py"


rule plot_concordance:
    input:
        "work_folder/analysis/assay_concordance/concordance_any_pos.tsv"
    output:
        "work_folder/analysis/assay_concordance/plots/concordance.png"
    log:
        "logs/analysis/assay_concordance/plots/concordance.log"
    script:
        "scripts/plot_concordance.R"
