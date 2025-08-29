


rule get_go_do_enrichmnets:
    input:
        hippie_degree="work_folder/degree/full_hippie.csv",
        summed_probability="work_folder/degree/flat_summed.csv",
        threshold_1="work_folder/degree/flat_min.1.csv",
        threshold_2="work_folder/degree/flat_min.2.csv"
    output:
        distribution_plot = "work_folder/plots/degree/Distribution.png",
        go_enrichment = "work_folder/plots/degree/GO_enrichment.png",
        do_enrichment = "work_folder/plots/degree/DO_enrichment.png"
    run:
