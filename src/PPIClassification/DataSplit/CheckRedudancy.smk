def read_fasta(fasta_filename):
    gene_name_seq_dict = dict()
    with open(fasta_filename) as f:
        gene_name = ""
        for line in f:
            if line[0] == ">":
                if gene_name:
                    gene_name = gene_name.groups()[0]
                    gene_name_seq_dict[gene_name] = seq
                gene_name = re.search(" GN=([A-Za-z/0-9-]+) ",line)
                seq = ""
            else:
                seq += line.strip()
    return gene_name_seq_dict

rule subset_fasta:
    # NOTE: fasta header lost
    input:
        partition=f"work_folder{pn}/subsets/{{subset}}/genes/genes_{{dataset}}_{{partition_name}}.txt",
        fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{dataset}}_{{partition_name}}.fasta"
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


# rule cdhit_2d:
#     params:
#         cdhit_location=config["cdhit_location"],
#         identity_threshold=0.5
#     input:
#         fasta_train=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}_{{partition_name}}.fasta"
#         fasta_test=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}_{{partition_name}}.fasta"
#     output:
        
#     shell:
#         """
#         {params.cdhit_location} -i {input.fasta_train} -i2 {input.fasta_} -o {output.pos} -o2 {output.neg} -c {params.identity_threshold} -n 5 -d 0 -M 16000 -T 8
#         """