rule get_cell_line_prey_plot:
    params:
        or_n = 50,
        script_location = "src/Plotting/CellLine/plot_prey_probability.R"
    input:
        cl_prey = "work_folder/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv"
    output:
        plot = "work_folder/plots/cell_line_prey.png"
    shell:
        """
        Rscript {params.script_location} \
             {input.cl_prey} \
             {params.or_n} \
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


