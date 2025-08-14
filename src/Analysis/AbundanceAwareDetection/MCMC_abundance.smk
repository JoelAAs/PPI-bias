import re
import pandas as pd
import glob
import os
import numpy as np
import math

def get_bait_parameters(wildcards):
    BAIT_FOLDER = checkpoints.batch_tests.get().output[0]
    bait_files = glob.glob(BAIT_FOLDER + "/*")
    batches = [
        re.search("/batch_(.*).csv$",f).groups()[0] for f in bait_files
    ]
    expected = [
        f"work_folder/analysis/abundance_aware/parameters_{wildcards.model}/batch_{current_batch}_parameters.csv"
        for current_batch in batches
    ]

    return expected


def get_cell_line_total(wildcards):
    BAIT_FOLDER = checkpoints.estimate_bait_interaction.get().output[0]
    bait_files = glob.glob(BAIT_FOLDER + "/*")

    return bait_files


rule get_shared_detections:
    input:
        all_genes="data/normalised_log_ra.csv"
    output:
        detectable_genes="work_folder/analysis/abundance/detectable_genes.csv"
    run:
        cl_abundances = pd.read_csv(input.all_genes,sep="\t")
        genes = cl_abundances.select_dtypes(np.float64).columns

        with open(output.detectable_genes,"w") as w:
            _ = [w.write(gene + "\n") for gene in genes]


checkpoint estimate_bait_interaction:
    params:
        min_tested=3,
        cellines=["CVCL_0030", "CVCL_0291", "CVCL_0063"]
    input:
        detectable_genes="work_folder/analysis/abundance/detectable_genes.csv",
        cl_specific_interactions="data/CL_annotated_bait_prey.csv"
    output:
        bait_folder=directory("work_folder/analysis/abundance_aware/baits")
    run:
        os.makedirs(output.bait_folder,exist_ok=True)
        PPI_interactions = pd.read_csv(input.cl_specific_interactions,sep="\t")
        with open(input.detectable_genes,"r") as f:
            detectable_genes = {l.strip() for l in f}

        PPI_interactions = PPI_interactions[
            (PPI_interactions["gene_name_bait"] != PPI_interactions["gene_name_prey"]) &
            (PPI_interactions["cl_id"].isin(params.cellines))
            ]
        PPI_interactions["study_id"] = PPI_interactions[
            ["pubmed_id", "detection_method"]
        ].apply(lambda row: "-".join(map(str,row)),axis=1)
        n_studies = PPI_interactions.groupby(["gene_name_bait", "cl_id"],as_index=False
        )["study_id"].nunique()
        n_studies = n_studies.rename(
            {
                "study_id": "n_tested"
            },axis=1
        )
        n_bait_prey_tests = PPI_interactions.groupby(
            ["gene_name_bait", "gene_name_prey", "cl_id"])["study_id"].nunique()

        total_studies = n_studies.groupby("gene_name_bait")["n_tested"].sum()
        baits_to_keep = total_studies[total_studies > params.min_tested].index.values
        n_studies = n_studies[n_studies["gene_name_bait"].isin(baits_to_keep)]

        for bait_name in n_studies["gene_name_bait"].unique():
            cl_ss = n_studies[n_studies["gene_name_bait"] == bait_name]
            with open(f"{output.bait_folder}/{bait_name}.csv","w") as w:
                w.write("\t".join(
                    [
                        "gene_name_bait",
                        "gene_name_prey",
                        "n_observed",
                        "n_tested",
                        "cl_id"
                    ]
                ) + "\n")
                for i, (_, cl_id, n_studies_bait) in cl_ss.iterrows():
                    for possible_prey in detectable_genes:
                        try:
                            interaction_observations = n_bait_prey_tests[bait_name][possible_prey][cl_id]
                        except KeyError:
                            interaction_observations = 0

                        w.write("\t".join(
                            map(str,[
                                bait_name,
                                possible_prey,
                                interaction_observations,
                                n_studies_bait,
                                cl_id
                            ])
                        ) + "\n")


rule get_shared_tests:
    input:
        bait_parameters=get_cell_line_total
    output:
        unique_prey_test_combinations="work_folder/analysis/abundance_aware/uniquely_tested_prey.csv"
    run:
        all_dfs = [
            pd.read_csv(bait_file,sep="\t") for bait_file in input.bait_parameters
        ]
        all_dfs = pd.concat(all_dfs)
        pivot_observed = all_dfs.pivot_table(index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',values='n_observed',fill_value=0)
        pivot_observed = pivot_observed.rename({
            c_name: f"n_observed_{c_name}" for c_name in pivot_observed
        },axis=1)

        pivot_tested = all_dfs.pivot_table(index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',values='n_tested',fill_value=0)
        pivot_tested = pivot_tested.rename({
            c_name: f"n_tested_{c_name}" for c_name in pivot_tested
        },axis=1)

        pivot_df = pivot_tested.join(pivot_observed)
        id_cols = pivot_df.columns.values.tolist()

        unique_tests = pivot_df.reset_index()
        del unique_tests["gene_name_bait"]
        unique_tests = unique_tests[~unique_tests[["gene_name_prey"] + id_cols].duplicated()]
        unique_tests.to_csv(output.unique_prey_test_combinations,sep="\t",index=False)


checkpoint batch_tests:
    params:
        batch_size=config["prey_n_MCMC_batch"]
    input:
        unique_prey_test_combinations="work_folder/analysis/abundance_aware/uniquely_tested_prey.csv"
    output:
        batch_folder=directory("work_folder/analysis/abundance_aware/batched_prey_tests")
    run:
        os.makedirs(output.batch_folder,exist_ok=True)
        unique_per_prey_df = pd.read_csv(input.unique_prey_test_combinations,sep="\t")
        i = 0
        batch = 0
        nrow = unique_per_prey_df.shape[0]
        while nrow > i:
            unique_per_prey_df.iloc[i:(i + params.batch_size)].to_csv(
                f"{output.batch_folder}/batch_{batch}.csv",sep="\t",index=False
            )
            i += params.batch_size
            batch += 1


rule fit_parameters_abundance:
    params:
        workers=config["MCMC_workers"],
        batch_size=config["batch_per_write"],# for now
        samples=config["mcmc_samples"],
        burin_samples=config["mcmc_burin"]
    input:
        bait="work_folder/analysis/abundance_aware/batched_prey_tests/batch_{batch}.csv",
        abundance_cell_lines="data/normalised_log_ra.csv"
    output:
        bait_parameters="work_folder/analysis/abundance_aware/parameters_abundance/batch_{batch}_parameters.csv"
    shell:
        """
        python src/Analysis/AbundanceAwareDetection/fit_model_mp.py \
            --prey_tested {input.bait} \
            --abundance_cell_lines {input.abundance_cell_lines} \
            --abundance 1 \
            --bait_output {output.bait_parameters} \
            --workers {params.workers} \
            --batch_size {params.batch_size} \
            --samples {params.samples} \
            --burin_samples {params.burin_samples}
        """


rule fit_parameters_pod:
    params:
        workers=config["MCMC_workers"],
        batch_size=config["batch_per_write"],
        samples=config["mcmc_samples"],
        burin_samples=config["mcmc_burin"]
    input:
        bait="work_folder/analysis/abundance_aware/batched_prey_tests/batch_{batch}.csv",
        abundance_cell_lines="data/normalised_log_ra.csv"
    output:
        bait_parameters="work_folder/analysis/abundance_aware/parameters_pod/batch_{batch}_parameters.csv"
    shell:
        """
        python src/Analysis/AbundanceAwareDetection/fit_model_mp.py \
            --prey_tested {input.bait} \
            --abundance_cell_lines {input.abundance_cell_lines} \
            --abundance 0 \
            --bait_output {output.bait_parameters} \
            --workers {params.workers} \
            --batch_size {params.batch_size} \
            --samples {params.samples} \
            --burin_samples {params.burin_samples}
        """


rule aggregate:
    params:
        selected_cell_lines=config["selected_cell_lines"]
    input:
        bait_parameters=get_bait_parameters
    output:
        aggregate_parameters="work_folder/analysis/abundance_aware/parameters_{model}/all_parameters.csv"
    run:
        betas = [f"beta_prediction_{c}_{suffix}" for c in
                 params.selected_cell_lines for suffix in ["mean", "sd"]]
        with open(output.aggregate_parameters,"w") as w:
            w.write("\t".join(
                ['gene_name_prey'] +
                [f"n_tested_{c}" for c in params.selected_cell_lines] +
                [f"n_observed_{c}" for c in params.selected_cell_lines] +
                betas +
                [
                    "beta_bait_mean",
                    "beta_bait_sd",
                    "beta_bait_low_ci",
                    "beta_bait_high_ci",
                    "n_divergences"
                ]
            ) + "\n")
            for bait_parameters in input.bait_parameters:
                with open(bait_parameters,"r") as f:
                    for line in f:
                        w.write(line)


rule get_bait_prey_pairs:
    input:
        baits_preys=get_cell_line_total,
        models="work_folder/analysis/abundance_aware/parameters_{model}/all_parameters.csv"
    output:
        all_bait_prey_models="work_folder/analysis/abundance_aware/bait_prey_{model}.csv"
    run:
        all_bait_dfs = [pd.read_csv(bp,sep="\t") for bp in input.baits_preys]
        all_bait_dfs = pd.concat(all_bait_dfs)
        pivot_observed = all_bait_dfs.pivot_table(index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',values='n_observed',fill_value=0)
        pivot_observed = pivot_observed.rename({
            c_name: f"n_observed_{c_name}" for c_name in pivot_observed
        },axis=1)

        pivot_tested = all_bait_dfs.pivot_table(index=['gene_name_bait', 'gene_name_prey'],
            columns='cl_id',values='n_tested',fill_value=0)
        pivot_tested = pivot_tested.rename({
            c_name: f"n_tested_{c_name}" for c_name in pivot_tested
        },axis=1)

        pivot_df = pivot_tested.join(pivot_observed)
        id_cols = pivot_df.columns.values.tolist()
        unique_tests = pivot_df.reset_index()
        unique_tests[id_cols] = unique_tests[id_cols].astype(int)

        model_params = pd.read_csv(input.models,sep="\t")
        full = unique_tests.merge(model_params,on=["gene_name_prey"] + id_cols)
        full.to_csv(output.all_bait_prey_models,sep="\t",index=False)


rule get_negatome_HCI:
    params:
        top_limit=0.5,
        bot_limit=0.05,
        selected_cell_lines=config["selected_cell_lines"]
    input:
        all_bait_prey_models="work_folder/analysis/abundance_aware/bait_prey_{model}.csv"
    output:
        negatome="work_folder/analysis/abundance_aware/negatome_{model}.csv",
        hci="work_folder/analysis/abundance_aware/HCI_{model}.csv"
    run:
        model_data = pd.read_csv(input.all_bait_prey_models,sep="\t")
        model_data["total_tested"] = model_data[[f"n_tested_{c}" for c in params.selected_cell_lines]].sum(axis=1)
        model_data["total_observed"] = model_data[[f"n_observed_{c}" for c in params.selected_cell_lines]].sum(axis=1)


        def mixture_mean(row, selected_cell_lines):
            mm = 0.0
            total_tested = row["total_tested"]
            for c in selected_cell_lines:
                value = row.get(f"beta_prediction_{c}_mean",np.nan)
                n_tests = row.get(f"n_tested_{c}",np.nan)
                if pd.notna(value) and pd.notna(n_tests) and total_tested != 0:
                    mm += float(value) * float(n_tests) / float(total_tested)
            return mm


        def mixture_var(row, samples):
            mv = 0.0
            c_mixture_mean = row["mixture_mean"]
            total_tested = row["total_tested"]

            for c in samples:
                mu = row.get(f"beta_prediction_{c}_mean",np.nan)
                sd = row.get(f"beta_prediction_{c}_sd",np.nan)
                n_tests = row.get(f"n_tested_CVCL_{c}",0)

                if pd.notna(mu) and pd.notna(sd) and n_tests != 0:
                    mv += (sd ** 2 + (mu - c_mixture_mean) ** 2) * n_tests / total_tested

            return mv


        model_data["mixture_mean"] = model_data.apply(lambda x: mixture_mean(x,params.selected_cell_lines),axis=1)
        model_data["mixture_var"] = model_data.apply(lambda x: mixture_var(x,params.selected_cell_lines),axis=1)

        model_data["lower_bound_pod"] = 2 ** model_data["mixture_mean"] * model_data["beta_bait_low_ci"]
        model_data["upper_bound_pod"] = 2 ** model_data["mixture_mean"] * model_data["beta_bait_high_ci"]

        negatome = model_data[model_data["upper_bound_pod"] < -math.log(1/params.bot_limit - 1)]
        negatome.to_csv(output.negatome, sep="\t", index=False)

        hci = model_data[model_data["lower_bound_pod"] > -math.log(1/params.top_limit - 1)]
        hci.to_csv(output.hci, sep="\t", index=False)