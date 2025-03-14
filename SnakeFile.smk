configfile: "config.yaml"

include: "src/FormatingFiltering.smk"

rule all:
    input:
        expand(
            "work_folder/intact/method_subset/interaction_counts/pair_count_{method}.csv",
            method = config["methods"]
        )
