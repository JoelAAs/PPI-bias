rule get_cell_line_prey_plot:
    params:
        part_or_difference_cutoff = config["part_or_difference_cutoff"],
        script_location = "src/Plotting/CellLine/plot_prey_probability.R"
    input:
        cl_prey = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv"
    output:
        plot = "work_folder/plots/cell_line_prey.png"
    shell:
        """
        Rscript {params.script_location} \
             {input.cl_prey} \
             {params.part_or_difference_cutoff} \
             {output.plot}
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
    shell:
        """
        Rscript {params.script_location} 
        """


rule plot_go:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_go_accumulation.R"
    input:
        greater_go="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_greater.csv",
        lesser_go="work_folder/analysis/GO/cumulative/POD_{data}_jaccard_lesser.csv"
    output:
        plot = "work_folder/plots/AccumulationPOD/go_{data}.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_go} {input.lesser_go} {wildcards.data} {output.plot}
        """


rule plot_colocalisation:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_colocalisation_accumulation.R"
    input:
        greater_colocalisation="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_greater.csv",
        lesser_colocalisation="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_lesser.csv",
    output:
        plot = "work_folder/plots/AccumulationPOD/colocalisation_{data}.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_colocalisation} {input.lesser_colocalisation} {wildcards.data} {output.plot}
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


