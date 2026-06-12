rule get_cell_line_prey_plot:
    params:
        part_or_difference_cutoff = config["part_or_difference_cutoff"],
        script_location = "src/Plotting/CellLine/plot_prey_probability.R"
    input:
        cl_prey = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv"
    output:
        plot = "work_folder/plots/cell_line_prey.png"
    log:
        "logs/plots/cell_line_prey.log"
    shell:
        """
        Rscript {params.script_location} \
             {input.cl_prey} \
             {params.part_or_difference_cutoff} \
             {output.plot} > {log} 2>&1
        """


rule get_localisation_y2h_ms_plot:
    """
    Localisation plot hardcoded for now
    """
    params:
        script_location = "src/Plotting/Method/plot_localisation_ms_y2h.R"
    input:
        ms_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
        y2h_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_y2h.csv"
    output:
        plot = "work_folder/plots/localisation_OR_y2h_ms.png"
    log:
        "logs/plots/localisation_OR_y2h_ms.log"
    shell:
        """
        Rscript {params.script_location} > {log} 2>&1
        """


rule plot_go:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_go_accumulation.R"
    input:
        greater_go="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_greater.csv",
        lesser_go="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_lesser.csv"
    output:
        plot_jaccard = "work_folder/plots/AccumulationPOD/go_{data}_jaccard.png",
        plot_accumulation= "work_folder/plots/AccumulationPOD/go_{data}_accumulation.png"
    log:
        "logs/plots/AccumulationPOD/go_{data}_jaccard.log"
    shell:
        """
        Rscript {params.script_location} {input.greater_go} {input.lesser_go} {wildcards.data} {output.plot_jaccard} {output.plot_accumulation} > {log} 2>&1
        """


rule plot_do:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_do_accumulation.R"
    input:
        greater_go="work_folder/analysis/DO/cumulative/POD_{data}_jaccard_greater.csv",
        lesser_go="work_folder/analysis/DO/cumulative/POD_{data}_jaccard_lesser.csv"
    output:
        plot_jaccard="work_folder/plots/AccumulationPOD/do_{data}_jaccard.png",
        plot_accumulation="work_folder/plots/AccumulationPOD/do_{data}_accumulation.png"
    log:
        "logs/plots/AccumulationPOD/do_{data}_jaccard.log"
    shell:
        """
        Rscript {params.script_location} {input.greater_go} {input.lesser_go} {wildcards.data} {output.plot_jaccard} {output.plot_accumulation} > {log} 2>&1
        """



rule plot_naive_colocalisation:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_colocalisation_accumulation.R"
    input:
        greater_colocalisation="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_greater.csv",
        lesser_colocalisation="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_lesser.csv",
    output:
        plot = "work_folder/plots/AccumulationPOD/colocalisation_{data}.png"
    log:
        "logs/plots/AccumulationPOD/colocalisation_{data}.log"
    shell:
        """
        Rscript {params.script_location} {input.greater_colocalisation} {input.lesser_colocalisation} {wildcards.data} {output.plot} > {log} 2>&1
        """

rule plot_matched_colocalisation:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_colocalisation_accumulation.R"
    input:
        matched_colocalisation_lesser="work_folder/analysis/localisation/study_match_probability/cumulative/POD_{data}_localisation_lesser.csv",
        matched_colocalisation_greater="work_folder/analysis/localisation/study_match_probability/cumulative/POD_{data}_localisation_greater.csv"
    output:
        plot = "work_folder/plots/AccumulationPOD/matched_colocalisation_{data}.png"
    shell:
        """
        Rscript {params.script_location} {input.matched_colocalisation_greater} {input.matched_colocalisation_lesser} {wildcards.data} {output.plot}
        """

rule plot_hydrophobicity:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_hydro_accumulation.R"
    input:
        greater_hydro="work_folder/analysis/hydrophobicity/cumulative/POD_{data}_netsurfp2_greater.csv",
        lesser_hydro="work_folder/analysis/hydrophobicity/cumulative/POD_{data}_netsurfp2_lesser.csv",
    output:
        plot = "work_folder/plots/AccumulationPOD/hydrophobicity_{data}.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_hydro} {input.lesser_hydro} {wildcards.data} {output.plot}
        """

rule plot_localisation_huri_bioplex:
    params:
        script_location = "src/Plotting/Method/plot_localisation_huri_bioplex.R"
    input:
        localisation_method = "work_folder/analysis/localisation/HuRI_vs_Bioplex.csv"
    output:
        png = "work_folder/plots/localisation/HuRI_bioplex.png"
    shell:
        """
        Rscript {params.script_location} {input.localisation_method} {output.png}
        """

rule plot_membrane_huri_bioplex:
    params:
        script_location = "src/Plotting/Membrane/plot_membrane_huri_bioplex.R"
    input:
        localisation_method = "work_folder/analysis/membrane/HuRI_vs_Bioplex_total_mean.csv"
    output:
        png = "work_folder/plots/membrane/HuRI_bioplex.png"
    shell:
        """
        Rscript {params.script_location} {input.localisation_method} {output.png}
        """

rule plot_mf_diff:
    params:
        script_location="src/Plotting/GO_terms/huri_vs_bioplex_mf_terms.R"
    input:
        compare = "work_folder/analysis/GO/huri_vs_bioplex_annotation.csv"
    output:
        png = "work_folder/plots/GO_terms/HuRI_bioplex_mf.png"
    conda:
        "do_enrichment"
    shell:
        """
        Rscript {params.script_location} {input.compare} {output.png}
        """
