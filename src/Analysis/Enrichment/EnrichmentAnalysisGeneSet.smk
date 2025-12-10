import pandas as pd
from IPython.core.pylabtools import retina_figure
from snakemake.io import expand
import numpy as np

rule get_enrichment:
    params:
        script = "src/Analysis/Enrichment/enrichment_degree.R"
    input:
        degree = f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        go_enrichment_bait = f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_bait_go.csv",
        go_enrichment_prey = f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_prey_go.csv",
        do_enrichment_bait = f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_bait_do.csv",
        do_enrichment_prey = f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_prey_do.csv"
    conda: "do_enrichment"
    shell:
        """
        Rscript {params.script} \
            {input.degree} \
            {output.go_enrichment_bait} \
            {output.go_enrichment_prey} \
            {output.do_enrichment_bait} \
            {output.do_enrichment_prey}
        """


def input_enrichments(wc, types, c_limits, c_ont):
    # should be ordered in config
    c_data = wc.data
    expected_input = []
    for c_type in types:
        if c_type == "HCNI":
            c_limit = c_limits[1]
        else:
            c_limit = c_limits[0]
        expected_input += expand(
            "work_folder{pn}/degree/enrichment/{data}_{type}_{limit}_{source}_{ont}.csv",
            pn=pn,
            data=c_data,
            type=c_type,
            limit=c_limit,
            source=["bait", "prey"],
            ont=c_ont
        )
    expected_input += expand(
        "work_folder{pn}/degree/enrichment/{data}_summed_{source}_{ont}.csv",
        pn=pn,
        data=c_data,
        source=["bait", "prey"],
        ont=c_ont
    )
    expected_input += expand(
        "work_folder{pn}/degree/enrichment/{data}_naive_{source}_{ont}.csv",
        pn=pn,
        data=c_data,
        source=["bait", "prey"],
        ont=c_ont
    )
    return expected_input


rule n_enriched_per_method:
    params:
        hci_limits = config["hci_limits"],
        hcni_tested = config["hcni_tested"]
    input:
        all_degree_enrichments = lambda wc: input_enrichments(
            wc, ["HCI", "delta", "HCNI"], [config["hci_limits"], config["hcni_tested"]], ["go","do"])
    output:
        n_enrichments = f"work_folder{pn}/degree/enrichment/significant_ontologies/{{data}}.csv"
    run:
        with open(output.n_enrichments, "w") as w:
            w.write("data\ttype\tsource\tlimit\tont\tn_enrichments\n")
            for c_enrichment in input.all_degree_enrichments:
                n_enrich = sum(1 for _ in open(c_enrichment, "r")) - 1
                base_name = c_enrichment.split("/")[-1]
                variables = base_name.split("_")
                ont = variables[-1].removesuffix(".csv")
                source = variables[-2]
                data = variables[0]
                type = variables[1]
                if type == "summed":
                    type = "HCI"
                    limit = "Expected"
                elif type=="naive":
                    type="HCI"
                    limit="None"
                else:
                    limit = variables[2]

                w.write(f"{data}\t{type}\t{source}\t{limit}\t{ont}\t{n_enrich}\n")


rule n_enriched_intact:
    input:
        intact_enrichments = expand(
            "work_folder{pn}/degree/enrichment/intact_{type}_{ont}.csv",
            pn=pn, type=["bait","prey"], ont = ["do", "go"]
        )
    output:
        n_enrichments=f"work_folder{pn}/degree/enrichment/intact_significant_ontologies/intact.csv"
    run:
        with open(output.n_enrichments,"w") as w:
            w.write("data\ttype\tsource\tlimit\tont\tn_enrichments\n")
            for c_enrichment in input.intact_enrichments:
                n_enrich = sum(1 for _ in open(c_enrichment,"r")) - 1
                base_name = c_enrichment.split("/")[-1]
                variables = base_name.split("_")
                ont = variables[-1].removesuffix(".csv")
                source = variables[-2]
                data = "Intact"
                type = "HCI"
                limit = "None"
                w.write(f"{data}\t{type}\t{source}\t{limit}\t{ont}\t{n_enrich}\n")


rule n_doids_gene_degree:
    input:
        degree = f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        doid_degree = f"work_folder{pn}/degree/doid/{{data_set_limit}}_doid.csv"
    params:
        script = "src/Analysis/Enrichment/get_ndoids.R"
    conda: "do_enrichment"
    shell:
        """
        Rscript {params.script} \
            {input.degree} \
            {output.doid_degree}
        """

def extreme_value_permutation_test(degree_str, input_doid_degree_file, permutation_naive, n_top_genes):
    naive_mean = permutation_naive.mean()
    df_degree = pd.read_csv(input_doid_degree_file,sep="\t")
    top_degrees = df_degree.nlargest(n_top_genes,f"degree_{degree_str}")
    mean_n_doids = top_degrees["n_doid"].mean()
    delta_mean = np.abs(mean_n_doids - naive_mean)

    n_extreme = np.sum(
        np.abs(permutation_naive-naive_mean) > delta_mean
    )
    p = n_extreme/len(permutation_naive)
    return mean_n_doids, p, naive_mean



rule test_top_degree_against_naive:
    params:
        permutations = 1000000, # probabiliy should make sure that same permutation isn't picked
        n_top_genes = 50,
        hci_limits = config["hci_limits"],
        hcni_limits = config["hcni_tested"]
    input:
        hci_degree = expand(
            "work_folder{pn}/degree/doid/{{data}}_HCI_{hci_limit}_doid.csv",
            pn=pn, hci_limit=config["hci_limits"]),
        hcni_degree= expand(
            "work_folder{pn}/degree/doid/{{data}}_HCNI_{hcni_limit}_doid.csv",
            pn=pn,hcni_limit=config["hcni_tested"]),

        naive_degree = f"work_folder{pn}/degree/doid/{{data}}_naive_doid.csv"
    output:
        doid_test = f"work_folder{pn}/degree/doid/{{data}}_tested.csv"
    run:
        df_naive_degree = pd.read_csv(input.naive_degree, sep="\t")
        top_naive_bait = df_naive_degree.nlargest(params.n_top_genes,f"degree_bait")
        top_naive_prey = df_naive_degree.nlargest(params.n_top_genes,f"degree_prey")
        naive_permute_dict = {
            "bait": np.array([top_naive_bait.sample(round(params.n_top_genes * 0.9))["n_doid"].mean()
            for _ in range(params.permutations)]),
            "prey": np.array([top_naive_prey.sample(round(params.n_top_genes * 0.9))["n_doid"].mean()
            for _ in range(params.permutations)]),
        }

        with open(output.doid_test,"w") as w:
            w.write(f"data\ttype\tlimit\tsource\thci_mean\tnaive_mean\tpermutation_p\n_permutations\n")
            for degree_type in ["bait", "prey"]:
                for hci_file, hci_limit in zip(input.hci_degree, params.hci_limits):
                    hci_mean, p_permuted, c_naive_mean =  extreme_value_permutation_test(
                        degree_type, hci_file, naive_permute_dict[degree_type], params.n_top_genes)
                    w.write(f"{wildcards.data}\tHCI\t{hci_limit}\t{degree_type}\t"
                            f"{hci_mean}\t{c_naive_mean}\t{p_permuted}\t{params.permutations}\n")

                for hcni_file, hcni_limit in zip(input.hcni_degree, params.hcni_limits):
                    hcni_mean, p_permuted, c_naive_mean =  extreme_value_permutation_test(
                        degree_type, hcni_file, naive_permute_dict[degree_type], params.n_top_genes)
                    w.write(f"{wildcards.data}\tHCNI\t{hcni_limit}\t{degree_type}\t"
                            f"{hcni_mean}\t{c_naive_mean}\t{p_permuted}\t{params.permutations}\n")




# rule get_bait_list:
#     # TODO: evaluate if this is used or useful
#     params:
#         other_ms_methods = [
#             "MI-0006",
#             "MI-0007",
#             "MI-0096",
#             "MI-0004",
#             "MI-0019"
#         ],
#         biotin_id = "MI-1314"
#     input:
#         intact = "data/bait_prey_publications.csv",
#         localisation_annotations = "data/gene_attribute_edges.txt",
#         uniprot_gene_name = "data/uniprot_to_gene_name.csv"
#     output:
#         bioid_baits="work_folder_{project}/enrichment_analysis/bait_lists/bioid_baits.csv",
#         ms_baits="work_folder_{project}/enrichment_analysis/bait_lists/ms_baits.csv",
#         shared_balanced="work_folder_{project}/enrichment_analysis/bait_lists/shared_baits.csv"
#     run:
#         intact_df = pd.read_csv(input.intact, sep="\t")
#
#         gene_name_to_uniprot = pd.read_csv(input.uniprot_gene_name,sep="\t")
#         intact_df = intact_df.merge(gene_name_to_uniprot, left_on="bait",right_on="uniprot_id")
#
#         bioid_ss = intact_df[intact_df["detection_method"] == params.biotin_id]
#         ms_ss = intact_df[intact_df["detection_method"].isin(params.other_ms_methods)]
#
#         bioid_ss["gene_name"].to_csv(output.bioid_baits, sep="\t", index=False)
#         ms_ss["gene_name"].to_csv(output.ms_baits,sep="\t",index=False)
#
#         bioid_bait_list = bioid_ss["gene_name"].tolist()
#         ms_bait_list = ms_ss["gene_name"].tolist()
#
#         shared_baits = set(bioid_bait_list) & set(ms_bait_list)
#         with open(output.shared_balanced, "w") as w:
#             w.write("gene_name\tbioid_data\n")
#             for bait_list, bioid_bool in zip([ms_bait_list, bioid_bait_list], [0,1]):
#                 for bait in bait_list:
#                     if bait in shared_baits:
#                         w.write(
#                             f"{bait}\t{bioid_bool}\n")
#
#
# rule bait_enrichment:
#     # TODO: evaluate id this is used or useful
#     params:
#         n_top_baits = 100
#     input:
#         bioid_baits = "work_folder_{project}/enrichment_analysis/bait_lists/bioid_baits.csv",
#         ms_baits = "work_folder_{project}/enrichment_analysis/bait_lists/ms_baits.csv"
#     output:
#         bioid_bait_enrichment_go_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_GO_bioid.csv",
#         bioid_bait_enrichment_do_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_DO_bioid.csv",
#         ms_bait_enrichment_go_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_GO_ms.csv",
#         ms_bait_enrichment_do_output = "work_folder_{project}/enrichment_analysis/enrichment/bait_enrichment_DO_ms.csv",
#         venn_plot_go = "work_folder_{project}/enrichment_analysis/plots/venn_diagram_goid.png",
#         venn_plot_doid = "work_folder_{project}/enrichment_analysis/plots/venn_diagram_doid.png",
#         venn_plot_bait= "work_folder_{project}/enrichment_analysis/plots/venn_diagram_bait.png"
#     shell:
#         """
#         Rscript src/enrichment_analysis.R \
#             {input.bioid_baits} \
#             {input.ms_baits} \
#             {output.bioid_bait_enrichment_go_output} \
#             {output.bioid_bait_enrichment_do_output} \
#             {output.ms_bait_enrichment_go_output} \
#             {output.ms_bait_enrichment_do_output} \
#             {output.venn_plot_go} \
#             {output.venn_plot_doid} \
#             {output.venn_plot_bait} \
#             {params.n_top_baits}
#         """
