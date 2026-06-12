import re
import time
import pandas as pd
import requests


def _batches(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


rule get_uniprot_sequences:
    input:
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
    log:
        "logs/protein_sequences/uniprot_canonical.log"
    run:
        df = pd.read_csv(input.gene_names, sep="\t")
        # Use base accession as query key; map base_ac -> entrez_id
        base_to_entrez = {
            row.uniprot_id.split("-")[0]: str(row.gene_name)
            for _, row in df.iterrows()
        }
        base_ids = list(base_to_entrez.keys())
        ac_pattern = re.compile(r'\|([A-Z0-9]+)\|')

        with open(output.fasta, "w") as w:
            for batch in _batches(base_ids, 200):
                query = " OR ".join(f"accession:{ac}" for ac in batch)
                resp = requests.get(
                    "https://rest.uniprot.org/uniprotkb/search",
                    params={"query": query, "format": "fasta", "size": len(batch) + 10}
                )
                if not resp.ok:
                    print(f"WARNING: batch request failed ({resp.status_code})")
                    continue

                for entry in resp.text.strip().split("\n>"):
                    if not entry:
                        continue
                    if not entry.startswith(">"):
                        entry = ">" + entry

                    header, *seq_lines = entry.split("\n")
                    ac_match = ac_pattern.search(header)
                    if not ac_match:
                        continue
                    ac = ac_match.group(1)
                    entrez = base_to_entrez.get(ac)
                    if entrez is None:
                        continue

                    # Replace GN= with the Entrez ID so embedding scripts key by it
                    if re.search(r' GN=\S+', header):
                        header = re.sub(r' GN=\S+', f' GN={entrez}', header)
                    else:
                        header = re.sub(r' PE=', f' GN={entrez} PE=', header)

                    w.write(header + "\n" + "\n".join(seq_lines) + "\n")

                time.sleep(0.05)


rule get_esmc_embeddings:
    params:
        model = config.get("esmc_model", "esmc_600m"),
        script_location = "src/PPIClassification/Embeddings/get_embeddings_esmc.py"
    input:
        fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
    output:
        embeddings = "work_folder/embeddings/canonical_ESMC.pt"
    log:
        "logs/embeddings/canonical_ESMC.log"
    container: "/beegfs/scratch/ieo7513/.snakemake/apptainer/huggingface-transformers-all-latest-gpu-latest.sif" # run with --apptainer-args="--nv"
    shell:
        """
        pip install --user --quiet esm
        python3 {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_output {output.embeddings} \
        > {log} 2>&1
        """


rule get_esm2_embeddings:
    params:
        model = config["embedding_model"],
        script_location = "src/PPIClassification/Embeddings/get_embeddings.py"
    input:
        fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
    output:
        embeddings = "work_folder/embeddings/canonical_ESM2.pt"
    log:
        "logs/embeddings/canonical_ESM2.log"
    container: "/beegfs/scratch/ieo7513/.snakemake/apptainer/huggingface-transformers-all-latest-gpu-latest.sif" # run with --apptainer-args="--nv"
    shell:
        """
        python3 {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_output {output.embeddings} \
        > {log} 2>&1
        """


rule get_mean_max_features:
    input:
        embeddings = "work_folder/embeddings/canonical_{esm_model}.pt"
    output:
        protein_embeddings="work_folder/embeddings/canonical_{esm_model}_mean_max.csv.gz",
    script:
        "get_mean_max_embeddings.py"
