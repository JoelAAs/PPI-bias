import datetime
import re
import pandas as pd
import requests

def get_sp_uniprot_gene_name(gene_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND taxonomy_id:9606 AND reviewed:true",
        "format": "fasta",
        "size":3
    }
    response = requests.get(url, params=params)
    if response.status_code == 429:
        raise RuntimeError("Leave room to breath!")
    if not response.ok:
        return f"> QGN={gene_name} FAILED\n"

    fasta = response.text

    if len(fasta) < 5:
         return f"> QGN={gene_name}\n"
    else:
        hits = fasta.split("\n>")  # Someone put > in their description

        hit_keep = hits[0]
        for i, hit in enumerate(hits):
            sp_match = re.search(r' GN=([^ ]+)',hit)
            if sp_match and sp_match.groups()[0] == gene_name:
                hit_keep = hit
                if i != 0:
                    hit_keep = ">" + hit_keep
                break

        lines = hit_keep.split("\n")
        lines[0] = lines[0] + f" QGN={gene_name}"

        return "\n".join(lines) + "\n" #NOTE: there will be some extra empty lines fix later


rule get_all_canonical_sequences:
    input:
        intact_uniprot_genes = f"work_folder{pn}/gene_names/uniprot_to_gene_name.csv"
    output:
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp.fasta"
    log:
        f"logs{pn}/protein_sequences/gene_name_sp.log"
    run:
        gene_name_uniprot_df = pd.read_csv(
            input.intact_uniprot_genes, sep="\t"
        )
        gene_names = set(gene_name_uniprot_df["gene_name"].unique())

        with open(output.fasta, "w") as w:
            start = datetime.datetime.now()
            current_percent = 0
            for i, gene_name in enumerate(gene_names):
                if i/len(gene_names) > current_percent/100:
                    current_percent += 1
                    current_time = datetime.datetime.now()
                    delta_time = current_time-start
                    remaining_time = delta_time.seconds*(100-current_percent)/current_percent
                    msg = f"{current_percent} % done. Time eclipsed {int(delta_time.seconds/60)} m, {round(delta_time.seconds % 60)} s."
                    msg += f" Est. remaining {int(remaining_time/60)} m, {round(remaining_time % 60)} s"
                    print(msg)

                w.write(
                    get_sp_uniprot_gene_name(gene_name)
                )

rule swissprot_gn_to_intact_gn:
    input:
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp.fasta"
    output:
        intact_to_sp = f"work_folder{pn}/gene_names/uniprot_to_sp.csv",
        fasta_dedup = f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    log:
        f"logs{pn}/gene_names/uniprot_to_sp.log"
    run:
        with open(output.intact_to_sp, "w") as w:
            w.write("sp_gene_name\tintact_gene_name\n")
            with open(input.fasta, "r") as f:
                for line in f:
                    if line[0] == ">":
                        line = line.strip()
                        sp_match = re.search(r' GN=([^ ]+)', line)
                        if sp_match: # drop non-matching
                            sp_match = sp_match.groups()[0]
                            intact_match = re.search(r' QGN=([^"]+)', line).groups()[0]
                            if ";" not in intact_match: # drop multi-genes
                                w.write(f"{sp_match}\t{intact_match}\n")

        gene_names = pd.read_csv(output.intact_to_sp, sep="\t")
        gene_names["same"] = gene_names["sp_gene_name"] == gene_names["intact_gene_name"]
        duplicated_gene_names = gene_names[gene_names["sp_gene_name"].duplicated(keep=False)]
        all_duplicated = duplicated_gene_names["sp_gene_name"].values
        sp_keep = duplicated_gene_names[duplicated_gene_names["same"]]["sp_gene_name"].values
        rest = duplicated_gene_names[~duplicated_gene_names["sp_gene_name"].isin(sp_keep)]
        print(f"dropped {len(set(rest["sp_gene_name"]))} genes without clear gene annotation")

        written_protiens = set()

        with open(output.fasta_dedup, "w") as w:
            with open(input.fasta, "r") as f:
                write_sequence = False
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line[0] == ">":
                        sp_match = re.search(r' GN=([^ ]+)',line)
                        intact_match = re.search(r' QGN=([^"]+)',line).groups()[0]
                        if sp_match:
                            sp_match = sp_match.groups()[0]
                            if (sp_match not in all_duplicated or intact_match in sp_keep) and sp_match not in written_protiens:
                                written_protiens.add(sp_match)
                                write_sequence = True
                            else: 
                                write_sequence = False
                        else:
                            write_sequence = False
                    if write_sequence:
                        w.write(line + "\n")




rule get_esmc_embeddings:
    params:
        model = config.get("esmc_model", "esmc_600m"),
        script_location = "src/PPIClassification/Embeddings/get_embeddings_esmc.py"
    input:
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        embeddings = f"work_folder{pn}/embeddings/canonical_ESMC.pt"
    log:
        f"logs{pn}/embeddings/canonical_ESMC.log"
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
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        embeddings = f"work_folder{pn}/embeddings/canonical_ESM2.pt"
    log:
        f"logs{pn}/embeddings/canonical_ESM2.log"
    container: "/beegfs/scratch/ieo7513/.snakemake/apptainer/huggingface-transformers-all-latest-gpu-latest.sif" # run with --apptainer-args="--nv"
    shell:
        """
        python3 {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_csv {output.embeddings} \
        > {log} 2>&1
        """


rule get_mean_max_features:
    input:
        embeddings = f"work_folder{pn}/embeddings/canonical_{{esm_model}}.pt"
    output:
        protein_embeddings=f"work_folder{pn}/embeddings/canonical_{{esm_model}}_mean_max.csv.gz",
    script:
        "get_mean_max_embeddings.py"
