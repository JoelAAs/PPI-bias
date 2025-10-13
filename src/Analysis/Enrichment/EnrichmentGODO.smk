rule get_go_do_enrichment:
    # TODO: rscript is hardcoded it will fail now with project name
    """
    Look at GO/DO enrichment of top degrees among different adjusted networks 
    """
    input:
        hippie_degree=f"work_folder/{pn}/degree/full_hippie.csv",
        summed_probability=f"work_folder/{pn}/degree/flat_summed.csv",
        threshold_1=f"work_folder/{pn}/degree/flat_min.1.csv",
        threshold_2=f"work_folder/{pn}/degree/flat_min.2.csv"
    output:
        distribution_plot=f"work_folder/{pn}/plots/degree/Distribution.png",
        go_enrichment=f"work_folder/{pn}/plots/degree/GO_enrichment.png",
        do_enrichment=f"work_folder/{pn}/plots/degree/DO_enrichment.png",
        top_delta_genes=f"work_folder/{pn}/plots/degree/genes_top_delta.png",
        doid_vs_degree=f"work_folder/{pn}/plots/degree/doid_vs_deg.png"
    conda: "do_enrichment"
    shell:
        """
        Rscript src/Analysis/Enrichment/enrichment_degree.R
        """

