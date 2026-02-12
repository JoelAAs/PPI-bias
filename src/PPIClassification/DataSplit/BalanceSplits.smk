rule maxflow_splits:
    params:
        script_location="src/PPIClassification/DataSplit/max_flow.py",
        min_max_flow=65
    log: "logs/maxflow/{partition_name}/maxflow_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    output:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/maxflow/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
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
    params:
        script_location="src/PPIClassification/DataSplit/ILP.py",
        accepted_missmatch=2
    threads: 20
    log: "logs/ilp/{partition_name}/ilp_{settype}_{dataset}_{neg_limit}_poslim_{pos_limit}.log"
    input:
        set_pos=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{pos_limit}}_{{partition_name}}_pos.csv",
        set_neg=f"work_folder{pn}/subsets/{{settype}}/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    output:
        balanced_pos=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_pos.csv",
        balanced_neg=f"work_folder{pn}/subsets/{{settype}}/balanced/{{dataset}}_limit_{{neg_limit}}_poslim_{{pos_limit}}_{{partition_name}}_neg.csv"
    shell:
        """(

        python3 {params.script_location} \
            --positive_data {input.set_max_flow_pos} \
            --negative_data {input.set_max_flow_neg} \
            --balanced_positive {output.balanced_pos} \
            --balanced_negative {output.balanced_neg} \
            --accepted_error {params.accepted_missmatch} \
            --threads {threads} \
            --subset {wildcards.settype}
        ) >{log} 2>&1"""
