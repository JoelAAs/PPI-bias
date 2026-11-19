def naive_set_settings(wildcards):
    side = "positive" if wildcards.edge_type == "pos" else "negative"
    return config["naive_split"][side][wildcards.exp_set]


rule download_hippie:
    params:
        url = config["naive_split"]["hippie_url"]
    output:
        hippie = "work_folder/data/hippie/HIPPIE-current.txt.gz"
    log:
        "logs/data/hippie/HIPPIE-current.log"
    shell:
        """
        wget {params.url} -O {output.hippie} > {log} 2>&1
        """


rule format_hippie_positive:
    input:
        hippie = "work_folder/data/hippie/HIPPIE-current.txt.gz"
    output:
        edges = "work_folder/naive_split/undirectional/full/pos/hippie.tsv"
    log:
        "logs/naive_split/undirectional/full/pos/hippie.log"
    run:
        hippie = pd.read_csv(input.hippie, sep="\t", usecols=["uniprot_accession_A", "uniprot_accession_B"])
        hippie = hippie.rename(columns={"uniprot_accession_A": "bait", "uniprot_accession_B": "prey"})
        hippie[["bait", "prey"]].to_csv(output.edges, sep="\t", index=False)


rule get_naive_experimental_edges:
    params:
        settings = naive_set_settings,
        bait_col = f"{config['id_pattern']}_bait",
        prey_col = f"{config['id_pattern']}_prey"
    input:
        table = lambda wildcards: naive_set_settings(wildcards)["path"].format(network_type=wildcards.network_type)
    output:
        edges = "work_folder/naive_split/{network_type}/full/{edge_type}/{exp_set}.tsv"
    resources:
        mem_gb=40
    log:
        "logs/naive_split/{network_type}/full/{edge_type}/{exp_set}.log"
    run:
        settings = params.settings
        value_cols = list(set(["n_observed", settings["column"]]))
        df = pd.read_parquet(input.table, columns=[params.bait_col, params.prey_col] + value_cols)
        df = df.rename(columns={params.bait_col: "bait", params.prey_col: "prey"})
        edges = filter_on_threshold(df, settings, settings["threshold"])
        edges[["bait", "prey"]].to_csv(output.edges, sep="\t", index=False)


rule naive_train_test_split:
    params:
        test_fraction = config["naive_split"]["test_fraction"],
        seed = lambda wildcards: int(wildcards.seed)
    input:
        pos = "work_folder/naive_split/{network_type}/full/pos/{pos_set}.tsv",
    output:
        train = "work_folder/naive_split/seed_{seed}_train/{network_type}/full/pos/{pos_set}.tsv",
        test = "work_folder/naive_split/seed_{seed}_test/{network_type}/full/pos/{pos_set}.tsv"
    log:
        "logs/naive_split/{network_type}/seed_{seed}/{pos_set}.log"
    script:
        "scripts/naive_split.py"


rule balance_naive_sets:
    input:
        pos = "work_folder/naive_split/{set}/{network_type}/full/pos/{pos_set}.tsv",
        neg = "work_folder/naive_split/{network_type}/full/neg/{neg_set}.tsv"
    output:
        balanced = "work_folder/naive_split/{set}/{network_type}/balanced/{pos_set}_{neg_set}.tsv"
    resources:
        mem_gb=40
    log:
        "logs/naive_split/{network_type}/{set}/balanced/{pos_set}_{neg_set}.log"
    script:
        "scripts/balance_naive_sets.py"
    

rule get_naive_no_negative:
    input:
        balanced = "work_folder/naive_split/{set}/{network_type}/balanced/{pos_set}_{neg_set}.tsv",
    output:
        balanced_no = "work_folder/naive_split/{set}/{network_type}/balanced/{pos_set}_{neg_set}-no.tsv"
    resources:
        mem_gb=40
    log:
        "logs/naive_split/{network_type}/{set}/balanced/{pos_set}_{neg_set}-no.log"
    script:
        "scripts/naive_no_negative.py"


rule generete_sets_for_ppi_split_pipeline:
    input:
        pos = "work_folder/naive_split/undirectional/full/pos/{pos_set}.tsv",
        neg = "work_folder/naive_split/undirectional/full/neg/{neg_set}.tsv"
    output:
        pos = "work_folder/naive_split/matched_for_ppi-split/pos/{pos_set}_{neg_set}.csv",
        neg = "work_folder/naive_split/matched_for_ppi-split/neg/{pos_set}_{neg_set}.csv"
    script:
        "scripts/prepare_for_ppi-split_pl.py"
