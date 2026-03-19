rule maxflow_splits:
    # Directional balancing
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
    log: "logs/maxflow/maxflow_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    resources:
        mem_gb=100
    input:
        set_pos=f"work_folder{pn}/subsets/{{dataset}}_directional_limit_{{pos_limit}}_pos.pq",
        set_neg=f"work_folder{pn}/subsets/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.pq"
    output:
        set_pos=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_neg.csv",
        balance_data=f"work_folder{pn}/subsets/maxflow/balance_data/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_pos.csv"
    shell:
        """(
        python3 {params.script_location} \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --max_flow_positive {output.set_pos} \
            --max_flow_negative {output.set_neg} \
            --balance_file {output.balance_data}
        ) >{log} 2>&1"""


rule balance_undirectional:
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_undirectional_limit_{{pos_limit}}_{{partition_name}}_pos.pq",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_undirectional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq"
    output:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/undirectionalbalanced/{{dataset}}_undirectional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/undirectionalbalanced/{{dataset}}_undirectional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    shell:
        """
        python3 src/PPIClassification/DataSplit/balance_undirectional.py \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --output_positive {output.balanced_pos} \
            --output_negative {output.balanced_neg} 
        """


rule generate_random_negative_set:
    input:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/{{selection}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/{{selection}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"

    output:
        set_random_neg=f"work_folder{pn}/subsets/{{settype}}/{{selection}}/{{dataset}}_{{network_type}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}-random_neg.csv"
    shell:
        """
        python3 src/PPIClassification/DataSplit/random_negative.py \
            --positive_data {input.balanced_pos} \
            --negative_data {input.balanced_neg} \
            --random_negative_data {output.set_random_neg} \
            --network_type {wildcards.network_type}
        """