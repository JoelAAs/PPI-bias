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
