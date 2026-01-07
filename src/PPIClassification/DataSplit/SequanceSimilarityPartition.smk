import pandas as pd
from ..Embeddings.get_embeddings import read_fasta

rule blast_sequence_similarity:
    params:
        n_threads = 20
    input:
        fasta = f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    output:
        similarity_tsv = f"work_folder{pn}/embeddings/all_vs_all.tsv"
    shell:
        """
        makeblastdb -dbtype prot -in {input.fasta}
        blastp -query gene_name_sp.fasta -db gene_name_sp.fasta \
        -outfmt "6 qseqid qtitle sseqid stitle evalue bitscore"  \
        -max_hsps 1 -num_threads {params.n_threads} -out all_vs_all.tsv
        ## Eval > 10 not reported
        """


rule get_METIS_adjacency_list:
    params:
        normalize=True
    input:
        similarity_tsv = f"work_folder{pn}/embeddings/all_vs_all.tsv",
        aa_seq_fasta= f"work_folder{pn}/embeddings/gene_name_sp.fasta"
    output:
        ""
    run:
        gene_seq_dict = read_fasta(input.aa_seq_fasta)
        ava_blast_df = pd.read_csv(input.similarity_tsv, header=None, sep="\t")
        ava_blast_df.columns = ["from", "to", "eval", "bitscore"]
        ava_blast_df["from_gene"] = ava_blast_df["from"].apply(lambda x: x.split("|")[2].replace("_HUMAN", ""))
        ava_blast_df["to_gene"] = ava_blast_df["to"].apply(lambda x: x.split("|")[2].replace("_HUMAN", ""))
        ava_blast_df["from_length"] = ava_blast_df["from_gene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df["to_length"] = ava_blast_df["to_gene"].apply(lambda x: len(gene_seq_dict.get(x,"")))




