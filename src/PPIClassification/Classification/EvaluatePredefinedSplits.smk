

rule get_hippie_extra_sequences:
    input:
        splits=expand(
            "results_matched/hippie_{neg_set}-{state}/{split}.csv",
            neg_set=["fdr", "pod"],
            state=["cand", "nonobs"],
            split=["train", "val", "test_balanced", "test_realistic"],
        ),
        canonical_fasta="work_folder/protein_sequences/uniprot_canonical.fasta",
    output:
        fasta="work_folder/protein_sequences/hippie_extra.fasta",
    log:
        "logs/protein_sequences/hippie_extra.log",
    run:
        split_proteins = set()
        for f in input.splits:
            df = pd.read_csv(f, dtype={"protein1": str, "protein2": str})
            split_proteins |= set(df["protein1"]) | set(df["protein2"])

        with open(input.canonical_fasta) as fh:
            canonical = {line.split("|")[1] for line in fh if line.startswith(">")}

        missing = sorted(split_proteins - canonical)
        ac_pattern = re.compile(r'\|([A-Z0-9]+)\|')
        written = 0
        with open(log[0], "w") as logf, open(output.fasta, "w") as w:
            logf.write(f"{len(split_proteins)} proteins in hippie splits, {len(missing)} not in canonical fasta\n")
            for batch in _batches(missing, 100):
                resp = requests.get(
                    "https://rest.uniprot.org/uniprotkb/search",
                    params={"query": " OR ".join(f"accession:{ac}" for ac in batch), "format": "fasta", "size": len(batch) + 10},
                )
                resp.raise_for_status()
                written += _parse_and_write(resp.text, {ac: ac for ac in batch}, ac_pattern, w)
                time.sleep(0.05)
            logf.write(f"{written}/{len(missing)} sequences written\n")
        if written != len(missing):
            raise ValueError(f"Only {written}/{len(missing)} hippie accessions returned by UniProt")


rule get_hippie_extra_esm2_embeddings:
    params:
        model = config["embedding_model"],
        script_location = "src/PPIClassification/Embeddings/get_embeddings.py"
    input:
        fasta = "work_folder/protein_sequences/hippie_extra.fasta"
    output:
        embeddings = "work_folder/embeddings/hippie_extra_ESM2.pt"
    log:
        "logs/embeddings/hippie_extra_ESM2.log"
    container: "/beegfs/scratch/ieo7513/.snakemake/apptainer/huggingface-transformers-all-latest-gpu-latest.sif" # run with --apptainer-args="--nv"
    shell:
        """
        python3 {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_output {output.embeddings} \
        > {log} 2>&1
        """


rule get_hippie_extra_prostT5_embeddings:
    params:
        model = config.get("prostt5_model", "Rostlab/ProstT5"),
        script_location = "src/PPIClassification/Embeddings/get_embeddings_prostt5.py"
    input:
        fasta = "work_folder/protein_sequences/hippie_extra.fasta"
    output:
        embeddings = "work_folder/embeddings/hippie_extra_ProstT5.pt"
    log:
        "logs/embeddings/hippie_extra_ProstT5.log"
    container: "/beegfs/scratch/ieo7513/.snakemake/apptainer/huggingface-transformers-all-latest-gpu-latest.sif" # run with --apptainer-args="--nv"
    shell:
        """
        python3 {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_output {output.embeddings} \
        > {log} 2>&1
        """


rule get_hippie_extra_mean_max_features:
    input:
        embeddings = "work_folder/embeddings/hippie_extra_{esm_model}.pt"
    output:
        protein_embeddings = "work_folder/embeddings/hippie_extra_{esm_model}_mean_max.csv.gz"
    script:
        "../Embeddings/get_mean_max_embeddings.py"


rule join_hippie_extra_embeddings:
    input:
        canonical = "work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
        extra = "work_folder/embeddings/hippie_extra_{esm_model}_mean_max.csv.gz"
    output:
        protein_embeddings = "work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz"
    shell:
        """
        (zcat {input.canonical}; zcat {input.extra} | tail -n +2) | gzip > {output.protein_embeddings}
        """


rule predef_xgboost:
    input:
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        protein_embeddings = "work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz"
    output:
        params="work_folder/classification/predef/xgboost/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_parameters.txt",
        saved_model="work_folder/classification/predef/xgboost/model/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
    log:
        "logs/classification/predef/xgboost/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_model.log",
    threads: 10
    resources:
        mem_gb=40,
    params:
        random_state=1234,
    script:
        "scripts/ppi_classify_xgboost_naive.py"


rule predef_model_metrics:
    input:
        saved_model="work_folder/classification/predef/xgboost/model/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}.joblib",
        train="results_matched/{pos_set}_{neg_set}-{state}/train.csv",
        test="results_matched/{pos_set}_{neg_set}-{test_state}/test_balanced.csv",
        protein_embeddings = "work_folder/embeddings/predef_{esm_model}_mean_max.csv.gz"
    output:
        metrics="work_folder/classification/predef/xgboost/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
        pred="work_folder/classification/predef/xgboost/predictions/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_pred.tsv.gz",
    log:
        "logs/classification/predef/xgboost/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.log",
    threads: 1
    resources:
        mem_gb=8,
    script:
        "../ModelEvaluation/scripts/evaluate_predefined_split.py"


rule predef_all_metrics:
    input:
        metrics=lambda wc: expand(
            "work_folder/classification/predef/xgboost/metrics/{network_type}_{pos_set}_{neg_set}_{state}_{esm_model}_test_{test_state}_metrics.tsv",
            network_type="undirectional",
            esm_model = ["ProstT5", "ESM2"],
            pos_set=["hippie", "fdr", "pod"],
            neg_set=["fdr", "pod"],
            state=["cand", "nonobs"],
            test_state=["cand", "nonobs"]
        ),
    output:
        all_metrics="work_folder/classification/predef/xgboost/all_metrics.tsv",
    shell:
        "awk 'FNR > 1 || NR == 1' {input.metrics} > {output.all_metrics}"
