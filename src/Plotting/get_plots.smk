rule get_cell_line_prey_plot:
    params:
        part_or_difference_cutoff = config["part_or_difference_cutoff"],
        script_location = "src/Plotting/CellLine/plot_prey_probability.R"
    input:
        cl_prey = f"work_folder{pn}/inferred_search_space/analysis/cell_line/bait_wise_prey_plotting.csv"
    output:
        plot = f"work_folder{pn}/plots/cell_line_prey.png"
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
        ms_diff_localisation=f"work_folder{pn}/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
        y2h_diff_localisation=f"work_folder{pn}/inferred_search_space/analysis/localisation/diff_localisation_y2h.csv"
    output:
        plot = f"work_folder{pn}/plots/localisation_OR_y2h_ms.png"
    shell:
        """
        Rscript {params.script_location} 
        """


rule plot_go:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_go_accumulation.R"
    input:
        greater_go=f"work_folder{pn}/analysis/GO/cumulative/POD_{{data}}_jaccard_greater.csv",
        lesser_go=f"work_folder{pn}/analysis/GO/cumulative/POD_{{data}}_jaccard_lesser.csv"
    output:
        plot_jaccard = f"work_folder{pn}/plots/AccumulationPOD/go_{{data}}_jaccard.png",
        plot_accumulation= f"work_folder{pn}/plots/AccumulationPOD/go_{{data}}_accumulation.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_go} {input.lesser_go} {wildcards.data} {output.plot_jaccard} {output.plot_accumulation}
        """


rule plot_do:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_do_accumulation.R"
    input:
        greater_go=f"work_folder{pn}/analysis/DO/cumulative/POD_{{data}}_jaccard_greater.csv",
        lesser_go=f"work_folder{pn}/analysis/DO/cumulative/POD_{{data}}_jaccard_lesser.csv"
    output:
        plot_jaccard=f"work_folder{pn}/plots/AccumulationPOD/do_{{data}}_jaccard.png",
        plot_accumulation=f"work_folder{pn}/plots/AccumulationPOD/do_{{data}}_accumulation.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_go} {input.lesser_go} {wildcards.data} {output.plot_jaccard} {output.plot_accumulation}
        """



rule plot_naive_colocalisation:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_colocalisation_accumulation.R"
    input:
        greater_colocalisation=f"work_folder{pn}/analysis/localisation/cumulative/POD_{{data}}_localisation_greater.csv",
        lesser_colocalisation=f"work_folder{pn}/analysis/localisation/cumulative/POD_{{data}}_localisation_lesser.csv",
    output:
        plot = f"work_folder{pn}/plots/AccumulationPOD/colocalisation_{{data}}.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_colocalisation} {input.lesser_colocalisation} {wildcards.data} {output.plot}
        """

rule plot_matched_colocalisation:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_colocalisation_accumulation.R"
    input:
        matched_colocalisation_lesser=f"work_folder{pn}/analysis/localisation/study_match_probability/cumulative/POD_{{data}}_localisation_lesser.csv",
        matched_colocalisation_greater=f"work_folder{pn}/analysis/localisation/study_match_probability/cumulative/POD_{{data}}_localisation_greater.csv"
    output:
        plot = f"work_folder{pn}/plots/AccumulationPOD/matched_colocalisation_{{data}}.png"
    shell:
        """
        Rscript {params.script_location} {input.matched_colocalisation_greater} {input.matched_colocalisation_lesser} {wildcards.data} {output.plot}
        """

rule plot_hydrophobicity:
    params:
        script_location = "src/Plotting/AccumulationPOD/plot_hydro_accumulation.R"
    input:
        greater_hydro=f"work_folder{pn}/analysis/hydrophobicity/cumulative/POD_{{data}}_netsurfp2_greater.csv",
        lesser_hydro=f"work_folder{pn}/analysis/hydrophobicity/cumulative/POD_{{data}}_netsurfp2_lesser.csv",
    output:
        plot = f"work_folder{pn}/plots/AccumulationPOD/hydrophobicity_{{data}}.png"
    shell:
        """
        Rscript {params.script_location} {input.greater_hydro} {input.lesser_hydro} {wildcards.data} {output.plot}
        """

rule plot_localisation_huri_bioplex:
    params:
        script_location = "src/Plotting/Method/plot_localisation_huri_bioplex.R"
    input:
        localisation_method = f"work_folder{pn}/analysis/localisation/HuRI_vs_Bioplex.csv"
    output:
        png = f"work_folder{pn}/plots/localisation/HuRI_bioplex.png"
    shell:
        """
        Rscript {params.script_location} {input.localisation_method} {output.png}
        """

rule plot_membrane_huri_bioplex:
    params:
        script_location = "src/Plotting/Membrane/plot_membrane_huri_bioplex.R"
    input:
        localisation_method = f"work_folder{pn}/analysis/membrane/HuRI_vs_Bioplex_total_mean.csv"
    output:
        png = f"work_folder{pn}/plots/membrane/HuRI_bioplex.png"
    shell:
        """
        Rscript {params.script_location} {input.localisation_method} {output.png}
        """

rule plot_mf_diff:
    params:
        script_location="src/Plotting/GO_terms/huri_vs_bioplex_mf_terms.R"
    input:
        compare = f"work_folder{pn}/analysis/GO/huri_vs_bioplex_annotation.csv"
    output:
        png = f"work_folder{pn}/plots/GO_terms/HuRI_bioplex_mf.png"
    shell:
        """
        Rscript {params.script_location} {input.compare} {output.png}
        """
