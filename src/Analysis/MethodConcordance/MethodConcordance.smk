rule match_ms_y2h_detection:
    params:
        col_a=f"{config['id_pattern']}_bait",
        col_b=f"{config['id_pattern']}_prey"
    input:
        pod_ms="work_folder/analysis/POD/{network_type}/POD_ms.pq",
        pod_y2h="work_folder/analysis/POD/{network_type}/POD_y2h.pq"
    output:
        detection="work_folder/analysis/method_concordance/{network_type}/ms_y2h_detection.pq"
    log:
        "logs/analysis/method_concordance/{network_type}/ms_y2h_detection.log"
    script:
        "scripts/match_ms_y2h_detection.py"


rule test_shared_negatives:
    input:
        detection="work_folder/analysis/method_concordance/{network_type}/ms_y2h_detection.pq"
    output:
        fit="work_folder/analysis/method_concordance/{network_type}/shared_negatives_fit.rds",
        summary="work_folder/analysis/method_concordance/{network_type}/shared_negatives_summary.txt"
    log:
        "logs/analysis/method_concordance/{network_type}/shared_negatives.log"
    script:
        "scripts/test_shared_negatives.R"


rule plot_detection_percentage:
    input:
        detection="work_folder/analysis/method_concordance/{network_type}/ms_y2h_detection.pq",
        fit="work_folder/analysis/method_concordance/{network_type}/shared_negatives_fit.rds"
    output:
        mixed_effects_dist="work_folder/analysis/method_concordance/{network_type}/mixed_effects_dist.png",
        hub_adjusted_concordance="work_folder/analysis/method_concordance/{network_type}/hub_adjusted_concordance.png"
    log:
        "logs/analysis/method_concordance/{network_type}/plot_detection_percentage.log"
    script:
        "scripts/plot_detection_percetage.R"
