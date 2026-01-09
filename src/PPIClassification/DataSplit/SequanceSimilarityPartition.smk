import pandas as pd
import re
from ..Embeddings.get_embeddings import read_fasta

rule blast_sequence_similarity:
    params:
        n_threads = 45
    input:
        fasta = f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        similarity_tsv = f"work_folder{pn}/protein_sequences/similarity/all_vs_all.tsv"
    shell:
        """
        makeblastdb -dbtype prot -in {input.fasta} -title "Gene Name SP DB"
        blastp -query {input.fasta} -db {input.fasta} \
        -outfmt "6 qseqid stitle evalue bitscore"  \
        -max_hsps 1 -num_threads {params.n_threads} -out all_vs_all.tsv
        ## Eval > 10 not reported
        """


rule get_METIS_adjacency_list:
    params:
        normalize=True
    input:
        similarity_tsv = f"work_folder{pn}/protein_sequences/similarity/all_vs_all.tsv",
        aa_seq_fasta= f"work_folder{pn}/protein_sequences/gene_name_sp_dedup.fasta"
    output:
        similarity_edge_list =  f"work_folder{pn}/protein_sequences/similarity/avg_bitscore_all.fasta",

    run:
        gene_seq_dict = read_fasta(input.aa_seq_fasta)
        #mean_length = round(sum([len(s) for s in gene_seq_dict.values()])/len(gene_seq_dict))
        ava_blast_df = pd.read_csv(input.similarity_tsv, header=None, sep="\t")
        ava_blast_df.columns = ["qseqid", "stitle", "evalue", "bitscore"]
        id_dict = {
            x.split(" ")[0].split("|")[2].replace("_HUMAN", ""): re.search(r' GN=([^ ]+)',x).group(1) for
            x in ava_blast_df["stitle"].unique()
        }

        ava_blast_df["qgene"] = ava_blast_df["qseqid"].apply(lambda x: id_dict[x.split("|")[2].replace("_HUMAN", "")])
        ava_blast_df["sgene"] = ava_blast_df["stitle"].apply(lambda x: id_dict[x.split(" ")[0].split("|")[2].replace("_HUMAN", "")])
        ava_blast_df["q_length"] = ava_blast_df["qgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df["s_length"] = ava_blast_df["sgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df = ava_blast_df[ava_blast_df["qgene"] != ava_blast_df["sgene"]] # no loops
        ava_blast_df["bitscore_p_residue"] = ava_blast_df.apply(
            lambda x: x["bitscore"]/min([x["q_length"], x["s_length"]]), axis=1)
        ava_blast_df["edge_id"] = ava_blast_df[["qgene", "sgene"]].apply(lambda x: "-".join(sorted(x)), axis = 1)
        ava_blast_df = ava_blast_df[~ava_blast_df["edge_id"].duplicated(keep="first")]

        ava_blast_df[["qgene", "sgene", "bitscore_p_residue"]].to_csv(output.similarity_edge_list)

        with open(ouput.)
        nodes = range(1000)
        for node in nodes:
            nV nE



