import pandas as pd

rule get_bait_list:
    params:
        other_ms_methods = [
            "MI-0006",
            "MI-0007",
            "MI-0096",
            "MI-0004",
            "MI-0019"
        ],
        biotin_id = "MI-1314"
    input:
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/gene_attribute_edges.txt",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv"
    output:
        bioid_baits="work_folder_{project}/enrichment_analysis/bait_lists/bioid_baits.csv",
        ms_baits="work_folder_{project}/enrichment_analysis/bait_lists/ms_baits.csv",
        shared_balanced="work_folder_{project}/enrichment_analysis/bait_lists/shared_baits.csv"
    run:
        intact_df = pd.read_csv(input.intact, sep="\t")

        gene_name_to_uniprot = pd.read_csv(input.uniprot_gene_name,sep="\t")
        intact_df = intact_df.merge(gene_name_to_uniprot, left_on="bait",right_on="uniprot_id")

        bioid_ss = intact_df[intact_df["detection_method"] == params.biotin_id]
        ms_ss = intact_df[intact_df["detection_method"].isin(params.other_ms_methods)]

        bioid_ss["gene_name"].to_csv(output.bioid_baits, sep="\t", index=False)
        ms_ss["gene_name"].to_csv(output.ms_baits,sep="\t",index=False)

        bioid_bait_list = bioid_ss["gene_name"].tolist()
        ms_bait_list = ms_ss["gene_name"].tolist()

        shared_baits = set(bioid_bait_list) & set(ms_bait_list)
        with open(output.shared_balanced, "w") as w:
            w.write("gene_name\tbioid_data\n")
            for bait_list, bioid_bool in zip([ms_bait_list, bioid_bait_list], [0,1]):
                for bait in bait_list:
                    if bait in shared_baits:
                        w.write(
                            f"{bait}\t{bioid_bool}\n")


rule bait_enrichment:
    params:
        n_top_baits = 100
    input:
        bioid_baits = "work_folder_{project}/enrichment_analysis/bait_lists/bioid_baits.csv",
        ms_baits = "work_folder_{project}/enrichment_analysis/bait_lists/ms_baits.csv"
    output:
        bioid_bait_enrichment_go_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_GO_bioid.csv",
        bioid_bait_enrichment_do_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_DO_bioid.csv",
        ms_bait_enrichment_go_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_GO_ms.csv",
        ms_bait_enrichment_do_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_DO_ms.csv",
        venn_plot_go = "work_folder_{project}/enrichment_analysis/plots/venn_diagram_goid.png",
        venn_plot_doid = "work_folder_{project}/enrichment_analysis/plots/venn_diagram_doid.png",
        venn_plot_bait= "work_folder_{project}/enrichment_analysis/plots/venn_diagram_bait.png"
    shell:
        """
        Rscript src/enrichment_analysis.R \
            {input.bioid_baits} \
            {input.ms_baits} \
            {output.bioid_bait_enrichment_go_output} \
            {output.bioid_bait_enrichment_do_output} \
            {output.ms_bait_enrichment_go_output} \
            {output.ms_bait_enrichment_do_output} \
            {output.venn_plot_go} \
            {output.venn_plot_doid} \
            {output.venn_plot_bait} \
            {params.n_top_baits} 
        """
