
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
        -outfmt "6 qseqid sseqid evalue bitscore"  \
        -max_hsps 1 --num_threads {params.n_threads} -out all_vs_all.tsv
        """


rule get_METIS_adjacency_list:
    input:
        similarity_tsv = f"work_folder{pn}/embeddings/all_vs_all.tsv"
    output:
        ""
    run:
        #TODO: 