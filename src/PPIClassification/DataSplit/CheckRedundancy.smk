rule subset_fasta:
    # NOTE: fasta header lost
    input:
        partition=f"work_folder{pn}/subsets/{{subset}}/genes/genes_{{selected_data}}.txt",
        fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}.fasta"
    run:
        with open(input.partition, "r") as f:
            for line in f:
                partitiongene = [l.strip() for l in f]

        gene_seq_dict = read_fasta(input.fasta)
        with open(output.fasta, "w") as w:
            for gene in partitiongene:
                seq = gene_seq_dict.get(gene,"")
                if seq:
                    w.write(f">{gene}\n{seq}\n")


rule cdhit:
    params:
        cdhit_location=config["cdhit_location"],
        identity_threshold=0.4
    threads: 20
    input:
        fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}.fasta"
    output:
        sim_reduced_fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/cdhit/{{selected_data}}.fasta"
    shell:
        """
        {params.cdhit_location} -i {input.fasta} -o {output.sim_reduced_fasta} -c {params.identity_threshold} -n 2 -T {threads}
        """

rule cdhit_to_gene_list:
    input:
        sim_reduced_fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/cdhit/{{selected_data}}.fasta"
    output:
        gene_list=f"work_folder{pn}/subsets/{{subset}}/genes/cdhit/genes_{{selected_data}}.txt"
    run:
        with open(input.sim_reduced_fasta, "r") as f:
            genes = [line.strip()[1:] for line in f if line.startswith(">")]
        with open(output.gene_list, "w") as w:
            for gene in genes:
                w.write(gene + "\n")