
rule check_protein_promiscuity:
    input:
        bait_prey_pod = "work_folder/analysis/POD/directional/POD_{dataset}.pq",
    params:
        id_pattern = config["id_pattern"],
    output:
        bait_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_bait_detectability.tsv",
        prey_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_prey_detectability.tsv",
        variance_components = "work_folder/analysis/ProteinPromiscuity/{dataset}_variance_components.tsv",
        model_jls = "work_folder/analysis/ProteinPromiscuity/{dataset}_model.jls",
    log:
        "logs/analysis/ProteinPromiscuity/{dataset}.log",
    conda:
        "julia"
    threads:
        20
    shell:
        """
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} julia -t 1 src/Analysis/ProteinPromiscuity/scripts/model_protein_promiscuity.jl \
        {input.bait_prey_pod} {params.id_pattern} \
        {output.bait_detectability} {output.prey_detectability} \
        {output.variance_components} {output.model_jls} \
        > {log} 2>&1
        """

rule plot_cross_method_detectability:
    input:
        ms_bait_detectability = "work_folder/analysis/ProteinPromiscuity/ms_bait_detectability.tsv",
        ms_prey_detectability = "work_folder/analysis/ProteinPromiscuity/ms_prey_detectability.tsv",
        y2h_bait_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_bait_detectability.tsv",
        y2h_prey_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_prey_detectability.tsv"
    output:
        ms_y2h_plot = "work_folder/analysis/ProteinPromiscuity/plot/cross_method_detectability.png",
        bait_prey_plot = "work_folder/analysis/ProteinPromiscuity/plot/bait_prey_detectability.png"
    threads:
        1
    script:
        "scripts/plot_cross_method_detectability.R"


rule check_ms_sticky_proteins:
    input:
        bait_prey_pod = "work_folder/analysis/POD/directional/POD_ms.pq",
        model_jls = "work_folder/analysis/ProteinPromiscuity/ms_model.jls",
    output:
        ms_sticky_proteins = "work_folder/analysis/ProteinPromiscuity/ms_observed_vs_expected_prey.tsv"
    threads:
        5
    conda:
        "julia"
    script:
        "scripts/check_ms_sticky_proteins.jl"