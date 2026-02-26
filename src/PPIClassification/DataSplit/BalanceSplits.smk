rule maxflow_splits:
    # Directional balancing
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
        min_max_flow=65
    log: "logs/maxflow/{partition_name}/maxflow_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_directional_limit_{{pos_limit}}_{{partition_name}}_pos.pq",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq"
    output:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    shell:
        """(
        python3 {params.script_location} \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --max_flow_positive {output.set_pos} \
            --max_flow_negative {output.set_neg} \
            --min_max_flow {params.min_max_flow} \
            --subset {wildcards.settype}
        ) >{log} 2>&1"""

rule ilp:
    # Directional balancing
    params:
        script_location="src/PPIClassification/DataSplit/ILP.py",
        accepted_missmatch=2
    threads: 20
    log: "logs/ilp/{partition_name}/ilp_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_directional_limit_{{pos_limit}}_{{partition_name}}_pos.pq",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.pq"
    output:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/ilp/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/ilp/{{dataset}}_directional_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    shell:
        """(
        python3 {params.script_location} \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --balanced_positive {output.balanced_pos} \
            --balanced_negative {output.balanced_neg} \
            --accepted_error {params.accepted_missmatch} \
            --threads {threads}         
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
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_{{partition_name}}_pos.pq"
    output:
        set_neg=f"work_folder{pn}/subsets/{{settype}}/randomnegative/{{dataset}}_{{network_type}}_limit_{{pos_limit}}_{{partition_name}}-random_neg.pq"
    shell:
        """
        python3 src/PPIClassification/DataSplit/random_negative.py \
            --positive_data {input.set_pos} \
            --negative_data {output.set_neg} \
            --networkt_type {wildcards.network_type}
        """