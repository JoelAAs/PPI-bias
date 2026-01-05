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
    return response.text


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

rule get_esm_embeddings:
    params:
        model = config["embedding_model"],
        script_location = "src/PPI-classification/embeddings/get_embeddings.py"
    input:
        fasta = f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    output:
        embeddings_csv = f"work_folder{pn}/embeddings/canonical_embedding.csv.gz"
    conda:
        "huggingface"
    shell:
        """
        python {params.script_location} {input.fasta} {params.model} {output.embeddings_csv}
        """
