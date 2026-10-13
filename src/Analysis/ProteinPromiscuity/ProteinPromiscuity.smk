import random

def select_row_wise_study_files(study_folder, methods, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = study_folder

    # storage.fs() is not auto-prefixed, so add the prefix and retrieve the file ourselves
    prefix = workflow.storage_settings.default_storage_prefix
    query = f"{prefix.rstrip('/')}/{filename}" if prefix else filename
    storage_file = storage.fs(query)
    storage_object = storage_file.flags["storage_object"]
    storage_object.local_path().parent.mkdir(parents=True, exist_ok=True)
    storage_object.retrieve_object()
    ppi_df = pd.read_csv(storage_file, sep="\t")
    ppi_df = ppi_df[ppi_df["detection_method"].isin(methods)]
    ppi_df = ppi_df[
        ppi_df[f"{id_pattern}_bait"] != ppi_df[f"{id_pattern}_prey"]
    ]

    if ppi_df.empty:
        raise ValueError(f"No studies for method: {methods}")
    ppi_df = ppi_df[~ppi_df[[f"{id_pattern}_bait", f"{id_pattern}_prey", "pubmed_id"]].duplicated(keep="first")] # Remove isoforms, if considering gene names
    if remove_single:
        counts = ppi_df.groupby(["pubmed_id"], as_index=False).size()
        keep_pids = counts[counts["size"] != 1]["pubmed_id"]
        ppi_df = ppi_df[ppi_df["pubmed_id"].isin(keep_pids)]
    expected = [
        f"{STUDY_FOLDER}/{pid}_{method}.csv"
        for pid, method in ppi_df[["pubmed_id", "detection_method"]].drop_duplicates().itertuples(index=False)
    ]
    return expected

def get_row_wise_input_files(methods, id_pattern, filename, remove_single=True):
    STUDY_FOLDER = checkpoints.infer_experimental_search_space.get(cell_line="_method").output[0]
    STUDY_FOLDER = "work_folder" + STUDY_FOLDER.split("work_folder")[1]
    return select_row_wise_study_files(STUDY_FOLDER, methods, id_pattern, filename, remove_single)

rule get_row_wise_detection:
    input:
        studies = lambda wc: get_row_wise_input_files(config[wc.dataset], config["id_pattern"], "work_folder/formated/bait_prey_publications.csv"),
    output:
        row_wise_data = "work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
    log:
        "logs/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.log",
    threads:
        1
    run:
        studies = list(input.studies)
        study_folder = next((f for f in studies if os.path.isdir(f)), None)
        if study_folder:  # Checkpoint left unresolved (spawned job), select the studies ourselves
            studies = select_row_wise_study_files(
                study_folder, config[wildcards.dataset], config["id_pattern"],
                "work_folder/formated/bait_prey_publications.csv"
            )
        with open(output.row_wise_data, "w") as w:
            w.write("\t".join(["bait", "prey", "experiment", "detection"]) + "\n")
            for study in studies:
                with open(study, "r") as f:
                    next(f) # skip header
                    for l in f:
                        values = l.strip().split("\t")
                        bait = values[0]; prey = values[1]
                        detection = values[3]; exp = f"{values[5]}_{values[4]}"
                        w.write("\t".join([bait, prey, exp, detection]) + "\n")


rule fit_protein_promiscuity_model:
    # currently set to top 1000
    input:
        row_wise_data = ancient("work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv"),
    output:
        model_jls = "work_folder/analysis/ProteinPromiscuity/{dataset}_model.jls",
    log:
        "logs/analysis/ProteinPromiscuity/{dataset}_fit.log",
    conda:
        "julia"
    threads:
        15
    shell:
        """
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} julia -t 1 src/Analysis/ProteinPromiscuity/scripts/fit_protein_promiscuity_model.jl \
        {input.row_wise_data} {output.model_jls} 1000 \
        > {log} 2>&1
        """

rule evaluate_protein_promiscuity_model:
    input:
        model_jls = "work_folder/analysis/ProteinPromiscuity/{dataset}_model.jls",
    output:
        bait_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_bait_detectability.tsv",
        prey_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_prey_detectability.tsv",
        experiment_detectability = "work_folder/analysis/ProteinPromiscuity/{dataset}_experiment_detectability.tsv",
        variance_components = "work_folder/analysis/ProteinPromiscuity/{dataset}_variance_components.tsv",
    log:
        "logs/analysis/ProteinPromiscuity/{dataset}_evaluate.log",
    conda:
        "julia"
    threads: 
        1
    shell:
        """
        julia -t 1 src/Analysis/ProteinPromiscuity/scripts/evaluate_protein_promiscuity_model.jl \
        {input.model_jls} \
        {output.bait_detectability} {output.prey_detectability} \
        {output.experiment_detectability} {output.variance_components} \
        > {log} 2>&1
        """

rule plot_cross_method_detectability:
    input:
        ms_bait_detectability = "work_folder/analysis/ProteinPromiscuity/ms_bait_detectability.tsv",
        ms_prey_detectability = "work_folder/analysis/ProteinPromiscuity/ms_prey_detectability.tsv",
        y2h_bait_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_bait_detectability.tsv",
        y2h_prey_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_prey_detectability.tsv"
    output:
        ms_y2h_plot = "work_folder/analysis/ProteinPromiscuity/plot/cross_method_detectability.png",
        bait_prey_plot = "work_folder/analysis/ProteinPromiscuity/plot/bait_prey_detectability.png"
    threads:
        1
    script:
        "scripts/plot_cross_method_detectability.R"


rule check_ms_sticky_proteins:
    params:
        n_protein_subsample = 1000,
    input:
        row_wise_data = ancient("work_folder/analysis/ProteinPromiscuity/row_data/ms_row_wise_detection.tsv"),
        model_jls = "work_folder/analysis/ProteinPromiscuity/ms_model.jls",
    output:
        ms_sticky_proteins = "work_folder/analysis/ProteinPromiscuity/ms_dispersion_prey.tsv"
    log:
        "logs/analysis/ProteinPromiscuity/check_ms_sticky_proteins.log",
    threads:
        5
    conda:
        "julia"
    shell:
        """
        julia -t 1 src/Analysis/ProteinPromiscuity/scripts/check_ms_sticky_proteins.jl \
        {input.row_wise_data} {input.model_jls} \
        {params.n_protein_subsample} {output.ms_sticky_proteins} \
        > {log} 2>&1
        """

rule plot_ms_sticky_proteins:
    input:
        ms_sticky_proteins = "work_folder/analysis/ProteinPromiscuity/ms_dispersion_prey.tsv"
    output:
        dispersion_plot = "work_folder/analysis/ProteinPromiscuity/plot/ms_sticky_proteins_dispersion.png"
    threads:
        1
    script:
        "scripts/plot_ms_sticky_proteins.R"


rule get_ad_db_auto_activators:
    input:
        row_wise_data = ancient("work_folder/analysis/ProteinPromiscuity/row_data/y2h_row_wise_detection.tsv"),
        model_jls = "work_folder/analysis/ProteinPromiscuity/y2h_model.jls",
    params:
        n_protein_subsample = 1000,
    output:
        prey_auto_activators = "work_folder/analysis/ProteinPromiscuity/y2h_prey_auto_activators.tsv",
        bait_auto_activators = "work_folder/analysis/ProteinPromiscuity/y2h_bait_auto_activators.tsv"
    log:
        "logs/analysis/ProteinPromiscuity/get_y2h_auto_activators.log",
    threads:
        1
    conda:
        "julia"
    shell:
        """
        julia -t 1 src/Analysis/ProteinPromiscuity/scripts/get_y2h_auto_activators.jl \
        {input.row_wise_data} {input.model_jls} {params.n_protein_subsample} \
        {output.prey_auto_activators} {output.bait_auto_activators} \
        > {log} 2>&1
        """

rule plot_sticky_and_detectable_proteins:
    input:
        ms_sticky_proteins = "work_folder/analysis/ProteinPromiscuity/ms_dispersion_prey.tsv",
        ms_prey_detectability = "work_folder/analysis/ProteinPromiscuity/ms_prey_detectability.tsv",
    output:
        sticky_detectable_plot = "work_folder/analysis/ProteinPromiscuity/plot/sticky_and_detectable_proteins.png"
    threads:
        1
    script:
        "scripts/plot_sticky_and_detectable_proteins.R"


rule plot_y2h_auto_activators:
    input:
        prey_auto_activators = "work_folder/analysis/ProteinPromiscuity/y2h_prey_auto_activators.tsv",
        bait_auto_activators = "work_folder/analysis/ProteinPromiscuity/y2h_bait_auto_activators.tsv",
        y2h_prey_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_prey_detectability.tsv",
        y2h_bait_detectability = "work_folder/analysis/ProteinPromiscuity/y2h_bait_detectability.tsv",
    output:
        auto_activator_plot = "work_folder/analysis/ProteinPromiscuity/plot/y2h_auto_activators_detectability.png"
    threads:
        1
    script:
        "scripts/plot_y2h_auto_activators.R"


rule plot_degree_distribution:
    input:
        ms_bait_detectability = "work_folder/analysis/ProteinPromiscuity/ms_bait_detectability.tsv",
        ms_prey_detectability = "work_folder/analysis/ProteinPromiscuity/ms_prey_detectability.tsv",
    output:
        degree_distribution_plot = "work_folder/analysis/ProteinPromiscuity/plot/ms_degree_distribution.png"
    threads:
        1
    script:
        "scripts/plot_degree_distribution.R"


rule check_log_normality:
    params:
        min_tested = 300,
    input:
        row_wise_data = ancient("work_folder/analysis/ProteinPromiscuity/row_data/ms_row_wise_detection.tsv")
    output:
        log_normality_plot = "work_folder/analysis/ProteinPromiscuity/plot/ms_log_normality.png",
        log_normality_tsv = "work_folder/analysis/ProteinPromiscuity/ms_log_normality.tsv"
    log:
        "logs/analysis/ProteinPromiscuity/check_log_normality.log",
    threads:
        1
    shell:
        """
        python src/Analysis/ProteinPromiscuity/scripts/check_log_normal.py \
        {input.row_wise_data} {output.log_normality_plot} {output.log_normality_tsv} \
        {params.min_tested} \
        > {log} 2>&1
        """


rule rank_partner_effect:
    input:
        row_wise_data="work_folder/analysis/ProteinPromiscuity/row_data/{dataset}_row_wise_detection.tsv",
    output:
        ranked="work_folder/analysis/ProteinPromiscuity/autoactivators/{dataset}_{role}_partner_effect.tsv",
    params:
        min_pairs=20,
        rho0=0.3,
        n_simulations=3000,
    run:
        role = wildcards.role
        partner = "prey" if role == "bait" else "bait"
        df = pd.read_csv(input.row_wise_data, sep="\t")
        pair = df.groupby([role, partner], as_index=False).agg(
            n=("detection", "size"), k=("detection", "sum")
        )
        pair = pair[pair["n"] >= 2].copy()     # n==1 pairs are uninformative by construction
        pair["det"] = pair["k"] > 0            # detected at least once
        pair["full"] = pair["k"] == pair["n"]  # detected on every test
        tot = pair.groupby(role).agg(
            n_pairs=("n", "size"), n_obs=("n", "sum"), n_pos=("k", "sum"),
            n_det=("det", "sum"), n_full=("full", "sum"),
        )
        tot["p"] = tot["n_pos"] / tot["n_obs"]
        pair = pair.join(tot["p"], on=role)
        e = pair["n"] * pair["p"]
        v = pair["n"] * pair["p"] * (1 - pair["p"])
        pair["chi"] = ((pair["k"] - e) ** 2 / v).where(v > 0, 0.0)
        tot["D"] = pair.groupby(role)["chi"].sum() / (tot["n_pairs"] - 1)
        tot["f"] = (tot["n_full"] / tot["n_det"]).where(tot["n_det"] > 0) # Fraction of full replication

        # p = 0 or 1: every pair is uniform by construction, so D carries no information
        tot.loc[(tot["p"] == 0) | (tot["p"] == 1), "D"] = np.nan
        tot = tot[tot["n_pairs"] >= params.min_pairs]

        rng = np.random.default_rng(0)

        def uniformity_p(g):
            n, p, r, B = g["n"].to_numpy(), g["p"].iloc[0], params.rho0, params.n_simulations
            sim = rng.binomial(n, rng.beta(p * (1 - r) / r, (1 - p) * (1 - r) / r, size=(B, n.size)))
            q = sim.sum(axis=1, keepdims=True) / n.sum()  # re-estimate p as the observed stat does
            chi = ((sim - n * q) ** 2 / np.maximum(n * q * (1 - q), 1e-12)).sum(axis=1)
            return (1 + (chi <= g["chi"].sum()).sum()) / (B + 1)

        tot["p_vs_specific"] = (
            pair[pair[role].isin(tot.index[tot["D"].notna()])].groupby(role).apply(uniformity_p)
        )

        detected = tot["n_pos"] > 0
        tot["score"] = np.nan
        tot.loc[detected, "score"] = (-tot.loc[detected, "D"]).rank(
            pct=True
        ) + tot.loc[detected, "p"].rank(pct=True)
        tot.sort_values("score", ascending=False, na_position="last").to_csv(
            output.ranked, sep="\t"
        )
