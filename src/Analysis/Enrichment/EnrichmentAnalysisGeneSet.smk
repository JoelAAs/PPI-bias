import random

import pandas as pd
from snakemake.io import expand
import numpy as np
from scipy.stats import fisher_exact
from src.Analysis.Enrichment.go_annotation_support import get_go_genes, get_go_frequency
from src.Analysis.ProteomeComparison.paxdb_vs_hela import p_values

rule get_enrichment:
    params:
        script="src/Analysis/Enrichment/enrichment_degree.R"
    input:
        degree=f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        go_enrichment_bait=f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_bait_go.csv",
        go_enrichment_prey=f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_prey_go.csv",
        do_enrichment_bait=f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_bait_do.csv",
        do_enrichment_prey=f"work_folder{pn}/degree/enrichment/{{data_set_limit}}_prey_do.csv"
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
        hci_limits=config["hci_limits"],
        hcni_tested=config["hcni_tested"]
    input:
        all_degree_enrichments=lambda wc: input_enrichments(
            wc,["HCI", "delta", "HCNI"],[config["hci_limits"], config["hcni_tested"]],["go", "do"])
    output:
        n_enrichments=f"work_folder{pn}/degree/enrichment/significant_ontologies/{{data}}.csv"
    run:
        with open(output.n_enrichments,"w") as w:
            w.write("data\ttype\tsource\tlimit\tont\tn_enrichments\n")
            for c_enrichment in input.all_degree_enrichments:
                n_enrich = sum(1 for _ in open(c_enrichment,"r")) - 1
                base_name = c_enrichment.split("/")[-1]
                variables = base_name.split("_")
                ont = variables[-1].removesuffix(".csv")
                source = variables[-2]
                data = variables[0]
                type = variables[1]
                if type == "summed":
                    type = "HCI"
                    limit = "Expected"
                elif type == "naive":
                    type = "HCI"
                    limit = "None"
                else:
                    limit = variables[2]

                w.write(f"{data}\t{type}\t{source}\t{limit}\t{ont}\t{n_enrich}\n")


rule n_enriched_intact:
    input:
        intact_enrichments=expand(
            "work_folder{pn}/degree/enrichment/intact_{type}_{ont}.csv",
            pn=pn,type=["bait", "prey"],ont=["do", "go"]
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
        degree=f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        doid_degree=f"work_folder{pn}/degree/doid/{{data_set_limit}}_doid.csv"
    params:
        script="src/Analysis/Enrichment/get_ndoids.R"
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
        np.abs(permutation_naive - naive_mean) > delta_mean
    )
    p = n_extreme / len(permutation_naive)
    return mean_n_doids, p, naive_mean


rule test_top_degree_against_naive:
    params:
        permutations=1000000,# probability should make sure that same permutation isn't picked
        n_top_genes=50,
        hci_limits=config["hci_limits"],
        hcni_limits=config["hcni_tested"]
    input:
        hci_degree=expand(
            "work_folder{pn}/degree/doid/{{data}}_HCI_{hci_limit}_doid.csv",
            pn=pn,hci_limit=config["hci_limits"]),
        hcni_degree=expand(
            "work_folder{pn}/degree/doid/{{data}}_HCNI_{hcni_limit}_doid.csv",
            pn=pn,hcni_limit=config["hcni_tested"]),

        naive_degree=f"work_folder{pn}/degree/doid/{{data}}_naive_doid.csv",
        summed_degree=f"work_folder{pn}/degree/doid/{{data}}_summed.csv"
    output:
        doid_test=f"work_folder{pn}/degree/doid/{{data}}_tested.csv"
    run:
        df_naive_degree = pd.read_csv(input.naive_degree,sep="\t")
        top_naive_bait = df_naive_degree.nlargest(params.n_top_genes,"degree_bait")
        top_naive_prey = df_naive_degree.nlargest(params.n_top_genes,"degree_prey")
        # TODO: permutation not bootstrap, also the distribution may not be symmetric then doesnt work
        naive_permute_dict = {
            "bait": np.array([top_naive_bait.sample(round(params.n_top_genes * 0.9))["n_doid"].mean()
                              for _ in range(params.permutations)]),
            "prey": np.array([top_naive_prey.sample(round(params.n_top_genes * 0.9))["n_doid"].mean()
                              for _ in range(params.permutations)]),
        }

        with open(output.doid_test,"w") as w:
            w.write(f"data\ttype\tlimit\tsource\tdoid_mean\tnaive_doid_mean\tpermutation_p\tn_permutations\n")
            for degree_type in ["bait", "prey"]:
                for hci_file, hci_limit in zip(input.hci_degree,params.hci_limits):
                    hci_mean, p_permuted, c_naive_mean = extreme_value_permutation_test(
                        degree_type,hci_file,naive_permute_dict[degree_type],params.n_top_genes)
                    w.write(f"{wildcards.data}\tHCI\t{hci_limit}\t{degree_type}\t"
                            f"{hci_mean}\t{c_naive_mean}\t{p_permuted}\t{params.permutations}\n")

                for hcni_file, hcni_limit in zip(input.hcni_degree,params.hcni_limits):
                    hcni_mean, p_permuted, c_naive_mean = extreme_value_permutation_test(
                        degree_type,hcni_file,naive_permute_dict[degree_type],params.n_top_genes)
                    w.write(f"{wildcards.data}\tHCNI\t{hcni_limit}\t{degree_type}\t"
                            f"{hcni_mean}\t{c_naive_mean}\t{p_permuted}\t{params.permutations}\n")

                summed_mean, p_permuted, c_naive_mean = extreme_value_permutation_test(
                    degree_type,input.summed_degree,naive_permute_dict[degree_type],params.n_top_genes)
                w.write(f"{wildcards.data}\tHCI\texpected\t{degree_type}\t"
                        f"{summed_mean}\t{c_naive_mean}\t{p_permuted}\t{params.permutations}\n")


rule top_degree_get_doids:
    params:
        n_top_genes=50,
        script="src/Analysis/Enrichment/get_doid_frequency.R"
    input:
        degree=f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        doid_freq=expand("work_folder{pn}/degree/doid/freq/{{data_set_limit}}_count_{source}.csv",
            pn=pn,source=["bait", "prey"]),
        doid_annotated=expand("work_folder{pn}/degree/doid/freq/{{data_set_limit}}_annotated_{source}.csv",
            pn=pn,source=["bait", "prey"])
    conda: "do_enrichment"
    shell:
        """
        Rscript {params.script} \
            {input.degree} \
            {params.n_top_genes} \
            degree_bait \
            {output.doid_freq[0]} \
            {output.doid_annotated[0]} 
            
        Rscript {params.script} \
            {input.degree} \
            {params.n_top_genes} \
            degree_prey \
            {output.doid_freq[1]} \
            {output.doid_annotated[1]} 
        """


rule get_go_annotation:
    params:
        n_top_genes=50
    input:
        degree=f"work_folder{pn}/degree/{{data_set_limit}}.csv"
    output:
        go_frequency=expand(
            "work_folder{pn}/degree/GO/{{data_set_limit}}_count_{source}.csv",
            pn=pn,source=["bait", "prey"]
        )
    run:
        degree_df = pd.read_csv(input.degree,sep="\t")
        for source, output_file in zip(["bait", "prey"],output.go_frequency):
            top = degree_df.nlargest(params.n_top_genes,f"degree_{source}")
            go_dict = get_go_genes(top["gene_name"])
            go_frequency_df = get_go_frequency(go_dict)
            go_frequency_df.to_csv(output_file,sep="\t",index=False)


def test_go_binomial(go_file, naive_df, gos, n_tested):
    go_df = pd.read_csv(go_file,sep="\t")
    go_df = go_df.set_index("go_terms")
    ors = []
    p_values = []
    all_naive_obs = []
    all_set_obs = []
    gos = list(gos)
    for selected_go_term in gos:
        if selected_go_term in go_df.index:
            go_obs = go_df.loc[selected_go_term]["go_occurrence"]
        else:
            go_obs = 0
        if selected_go_term in naive_df.index:
            naive_obs = naive_df.loc[selected_go_term]["go_occurrence"]
        else:
            naive_obs = 0

        c_table = np.array(
            [[go_obs, n_tested - go_obs],
             [naive_obs, n_tested - naive_obs]]
        )
        odds_ratio, p = fisher_exact(c_table,alternative='two-sided')
        ors.append(odds_ratio)
        p_values.append(p)
        all_naive_obs.append(naive_obs)
        all_set_obs.append(go_obs)

    return [gos, ors, p_values, all_naive_obs, all_set_obs]

rule test_go_terms:
    params:
        n_tested_genes=50,
        min_observed=0.1,
        hci_limit=config["hci_limits"],
        hcni_tested=config["hcni_tested"]
    input:
        naive_go_frequency_df=f"work_folder{pn}/degree/GO/{{data}}_naive_count_{{source}}.csv.gz",
        hci_observations=expand(
            f"work_folder{pn}/degree/GO/{{data}}_HCI_{limit}_count_{{source}}.csv",
            pn=pn,limit=config["hci_limit"]),
        hcni_observations=expand(
            f"work_folder{pn}/degree/GO/{{data}}_HCNI_{limit}_count_{{source}}.csv",
            pn=pn,limit=config["hcni_tested"])
    output:
        test_csv = f"work_folder{pn}/degree/GO/tested/{{data}}_{{source}}.csv"
    run:
        naive_go_frequency_df = pd.read_csv(
            input.naive_go_frequency_df,sep="\t"
        )
        naive_go_frequency_df = naive_go_frequency_df.set_index("go_term")
        go_to_keep = set(naive_go_frequency_df[naive_go_frequency_df["go_frequency"] > params.min_observed]["go_term"])
        for alt_obs in input.hci_observations + input.hcni_observations:
            alt_obs_df = pd.read_csv(alt_obs,sep="\t")
            go_to_keep |= set(alt_obs_df[alt_obs_df["go_frequency"] > params.min_observed]["go_term"])

        go_set_zips = list(zip(
            input.hci_observations,params.hci_limit,["HCI"] * len(params.hci_limit)
        )) +  list(zip(
            input.hcni_observations,params.hcni_tested,["HCNI"] * len(params.hcni_tested)
        ))
        data_cols = ["go_term", "or", "p_value", "naive_obs", "set_obs"]
        id_cols = ["data","type" "limit", "source"]
        with open(output.test_csv, "w") as w:
            w.write("\t".join(id_cols + data_cols))

        for go_terms_file, limit, type_set in go_set_zips:
            go_results = test_go_binomial(go_terms_file, naive_go_frequency_df, go_to_keep, params.n_tested_genes)

            go_df = pd.DataFrame(go_results, columns=data_cols)
            go_df[id_cols] = [
                [f"{wildcards.data}"]*len(go_to_keep),
                [f"{type_set}"] * len(go_to_keep),
                [f"{limit}"] * len(go_to_keep),
                [f"{wildcards.source}"] * len(go_to_keep)
            ]
            go_df[id_cols + data_cols].to_csv(output.test_csv, sep="\t", mode="a", header=False, index=False)
