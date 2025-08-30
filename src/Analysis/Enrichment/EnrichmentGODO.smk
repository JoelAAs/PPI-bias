rule get_go_do_enrichment:
    input:
        hippie_degree="work_folder/degree/full_hippie.csv",
        summed_probability="work_folder/degree/flat_summed.csv",
        threshold_1="work_folder/degree/flat_min.1.csv",
        threshold_2="work_folder/degree/flat_min.2.csv"
    output:
        distribution_plot="work_folder/plots/degree/Distribution.png",
        go_enrichment="work_folder/plots/degree/GO_enrichment.png",
        do_enrichment="work_folder/plots/degree/DO_enrichment.png",
        top_delta_genes="work_folder/plots/degree/genes_top_delta.png",
        doid_vs_degree="work_folder/plots/degree/doid_vs_deg.png"
    shell:
        """
        Rscript src/Analysis/Enrichment/enrichment_degree.R
        """

