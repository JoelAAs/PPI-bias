import re

import pandas as pd
import requests
from urllib3.exceptions import RequestError
import time

def get_sp_uniprot_gene_name(gene_name):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND taxonomy_id:9606 AND reviewed:true",
        "format": "fasta",
        "size":1
    }
    response = requests.get(url, params=params)

    if not response.ok:
        print(response.status_code)
        raise RequestError(f"{gene_name} failed")

    fasta = response.text
    fasta_lines = fasta.split("\n")
    fasta_lines[0] += f"QGN: {gene_name}"
    return "\n".join(fasta_lines)


rule get_all_canonical_sequences:
    input:
        intact = f"work_folder{pn}/formated/bait_prey_publications.csv"
    output:
        fasta = f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    run:
        df_intact = pd.read_csv(
            input.intact, sep="\t"
        )
        gene_names = set(df_intact["gene_name_prey"]) | set(df_intact["gene_name_prey"])
        with open(output.fasta, "w") as w:
            for i, gene_name in enumerate(gene_names):
                print(f"{i}/{len(gene_names)} done.")
                time.sleep(0.1)
                w.write(
                    get_sp_uniprot_gene_name(gene_name)
                )

rule swissprot_gn_to_intact_gn:
    #TODO: This needs to be done in the beginning of intact. in Formating!
    input:
        fasta = f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    output:
        intact_to_sp = f"work_folder{pn}/embeddings/gn_sp2intact.csv",
        fasta_dedup = f"work_folder{pn}/embeddings/gene_name_sp_dedub.fasta"
    run:
        with open(outut.intact_to_sp, "w") as w:
            w.write("sp_gene_name\tintact_gene_name\n")
            with open(input.fasta, "r") as f:
                for line in f:
                    if line[0] == ">":
                        line = line.strip()
                        sp_match = re.search(r'GN=([^ ]+)', line).groups()[0]
                        intact_match = re.search(r'QGN: ([^"]+)', line).groups()[0]
                        w.write(f"{sp_match}\t{intact_match}\n")

        gene_names = pd.read_csv(output.intact_to_sp, sep="\t")
        gene_names["same"] = gene_names["sp_gene_name"] == gene_names["intact_gene_name"]
        duplicated_gene_names = gene_names[gene_names["sp_gene_name"].duplicated(keep=False)]
        all_duplicated = duplicated_gene_names["sp_gene_name"].values
        sp_keep = duplicated_gene_names[duplicated_gene_names["same"]]["sp_gene_name"].values
        rest = duplicated_gene_names[~duplicated_gene_names["sp_gene_name"].isin(sp_keep)]
        rest = rest[rest["sp_gene_name"].duplicated()]["intact_gene_name"].values

        with open(output.fasta_dedup, "w") as w:
            with open(input.fasta, "r") as f:
                write_sequence = False
                for line in f:
                    line = line.strip()
                    if line[0] == ">":
                        sp_match = re.search(r'GN=([^ ]+)',line).groups()[0]
                        intact_match = re.search(r'QGN: ([^"]+)',line).groups()[0]
                        if sp_match not in all_duplicated or intact_match in sp_keep or intact_match in rest:
                            write_sequence = True
                        else:
                            write_sequence = False
                    if write_sequence:
                        w.write(line + "\n")




rule get_esm_embeddings:
    params:
        model = config["embedding_model"],
        script_location = "src/PPIClassification/Embeddings/get_embeddings.py"
    input:
        fasta = f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    output:
        embeddings_csv = f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    conda:
        "huggingface"
    shell:
        """
        python {params.script_location} \
        --protein_fasta {input.fasta} \
        --model_name {params.model} \
        --embedding_csv {output.embeddings_csv}
        """
