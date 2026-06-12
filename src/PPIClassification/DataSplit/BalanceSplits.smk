
rule maxflow_splits:
    input:
        set_pos="work_folder/subsets/{dataset}_directional_full_{pos_limit}_pos.pq",
        set_neg="work_folder/subsets/{dataset}_directional_full_{neg_limit}_neg.pq",
    output:
        set_pos="work_folder/subsets/maxflow/{dataset}_directional_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
        set_neg="work_folder/subsets/maxflow/{dataset}_directional_limit_{neg_limit}_poslim_{pos_limit}_neg.csv",
        balance_data="work_folder/subsets/maxflow/balance_data/{dataset}_directional_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
    threads: 10
    resources:
        mem_gb=100,
    log:
        "logs/subsets/maxflow/{dataset}_directional_limit_{neg_limit}_poslim_{pos_limit}.log"
    # Directional balancing
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
    shell:
        """
        python3 {params.script_location} \
            --positive_data {input.set_pos} \
            --negative_data {input.set_neg} \
            --max_flow_positive {output.set_pos} \
            --max_flow_negative {output.set_neg} \
            --balance_file {output.balance_data} \
            --threads {threads} > {log} 2>&1
        """


rule balance_undirectional:
    input:
        set_pos="work_folder/subsets/{settype}/{dataset}_undirectional_limit_{pos_limit}_{partition_name}_pos.pq",
        set_neg="work_folder/subsets/{settype}/{dataset}_undirectional_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.pq",
    output:
        balanced_pos="work_folder/subsets/{settype}/undirectionalbalanced/{dataset}_undirectional_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_pos.csv",
        balanced_neg="work_folder/subsets/{settype}/undirectionalbalanced/{dataset}_undirectional_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}_neg.csv",
    log:
        "logs/subsets/{settype}/undirectionalbalanced/{dataset}_undirectional_limit_{neg_limit}_poslim_{pos_limit}_{partition_name}.log"
    script:
        "scripts/balance_undirectional.py"


rule generate_random_negative_set:
    input:
        balanced_pos="work_folder/subsets/{settype}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_pos.csv",
    output:
        set_random_neg="work_folder/subsets/{settype}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}-random_neg.csv",
    log:
        "logs/subsets/{settype}/{dataset}_{network_type}_limit_{neg_limit}_poslim_{pos_limit}_random_neg.log"
    script:
        "random_negative.py"