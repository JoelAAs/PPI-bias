import re
import time
import pandas as pd
import requests


def _batches(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def _parse_and_write(text, base_to_entrez, ac_pattern, w):
    written = 0
    for entry in text.strip().split("\n>"):
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

        if re.search(r' GN=\S+', header):
            header = re.sub(r' GN=\S+', f' GN={entrez}', header)
        else:
            header = re.sub(r' PE=', f' GN={entrez} PE=', header)

        w.write(header + "\n" + "\n".join(seq_lines) + "\n")
        written += 1
    return written


rule get_uniprot_sequences:
    params:
        id_pattern = config["id_pattern"]
    input:
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
    log:
        "logs/protein_sequences/uniprot_canonical.log"
    run:
        df = pd.read_csv(input.gene_names, sep="\t")
        base_to_entrez = {
            row.uniprot_id.split("-")[0]: row.uniprot_id.split("-")[0] if params.id_pattern == "uniprot_id" else str(row.gene_name)
            for _, row in df.iterrows()
        }
        ac_pattern = re.compile(r'\|([A-Z0-9]+)\|')

        with open(log[0], "w") as logf, open(output.fasta, "w") as w:
            # Bulk download reviewed entries (covers ~87% of accessions)
            logf.write(f"Downloading reviewed human proteome from UniProt...\n")
            logf.flush()
            resp = requests.get(
                "https://rest.uniprot.org/uniprotkb/stream",
                params={"query": "(reviewed:true) AND (organism_id:9606)", "format": "fasta"},
                stream=True,
            )
            if not resp.ok:
                raise RuntimeError(f"UniProt bulk download failed: HTTP {resp.status_code}")

            written = _parse_and_write(resp.text, base_to_entrez, ac_pattern, w)
            logf.write(f"Reviewed proteome: {written} sequences written\n")
            logf.flush()

            # Fetch remaining unreviewed accessions in batches of 100
            covered = set()
            for entry in resp.text.strip().split("\n>"):
                m = ac_pattern.search(entry)
                if m and m.group(1) in base_to_entrez:
                    covered.add(m.group(1))

            missing = [ac for ac in base_to_entrez if ac not in covered]
            logf.write(f"{len(missing)} unreviewed accessions to fetch in batches...\n")
            logf.flush()

            n_batches = (len(missing) + 99) // 100
            for i, batch in enumerate(_batches(missing, 100)):
                query = " OR ".join(f"accession:{ac}" for ac in batch)
                resp = requests.get(
                    "https://rest.uniprot.org/uniprotkb/search",
                    params={"query": query, "format": "fasta", "size": len(batch) + 10},
                )
                if not resp.ok:
                    logf.write(f"WARNING: batch {i+1}/{n_batches} failed (HTTP {resp.status_code})\n")
                    logf.flush()
                    continue
                batch_written = _parse_and_write(resp.text, base_to_entrez, ac_pattern, w)
                written += batch_written
                logf.write(f"Batch {i+1}/{n_batches}: {batch_written}/{len(batch)} written (total: {written})\n")
                logf.flush()
                time.sleep(0.05)

            logf.write(f"Done: {written}/{len(base_to_entrez)} sequences written to {output.fasta}\n")
            logf.flush()


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


# will only run with accenssions for now
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
