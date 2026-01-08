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
        hits = re.split("[(?:^|\n)]>", fasta)  # Someone put > in their description
        hit_keep = hits[0]
        for hit in hits:
            sp_match = re.search(r' GN=([^ ]+)',hit)
            if sp_match and sp_match.groups()[0] == gene_name:
                hit_keep = hit
                break

        lines = hit_keep.split("\n")
        lines[0] = lines[0] + f" QGN={gene_name}"

        return "\n".join(lines)


rule get_all_canonical_sequences:
    input:
        intact_uniprot_genes = f"work_folder{pn}/gene_names/uniprot_to_gene_name.csv"
    output:
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp.fasta"
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

        with open(output.fasta_dedup, "w") as w:
            with open(input.fasta, "r") as f:
                write_sequence = False
                for line in f:
                    line = line.strip()
                    if line[0] == ">":
                        sp_match = re.search(r' GN=([^ ]+)',line)
                        intact_match = re.search(r' QGN=([^"]+)',line).groups()[0]
                        if sp_match:
                            sp_match = sp_match.groups()[0]
                            if sp_match not in all_duplicated or intact_match in sp_keep:
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
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp_dedup.fasta"
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
