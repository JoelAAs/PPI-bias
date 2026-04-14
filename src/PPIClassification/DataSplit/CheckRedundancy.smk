rule subset_fasta:
    # NOTE: fasta header lost
    input:
        partition=f"work_folder{pn}/subsets/{{subset}}/genes/genes_{{selected_data}}.txt",
        fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}.fasta"
    log:
        f"logs{pn}/subsets/{{subset}}/genes/fasta/{{selected_data}}.log"
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
    log:
        f"logs{pn}/subsets/{{subset}}/genes/fasta/cdhit/{{selected_data}}.log"
    shell:
        """
        {params.cdhit_location}/cdhit -i {input.fasta} -o {output.sim_reduced_fasta} -c {params.identity_threshold} -n 2 -T {threads} > {log} 2>&1
        """

rule cdhit_redudance_between_subsets:
    params:
        cdhit_location=config["cdhit_location"],
        identity_threshold=0.4
    threads: 20
    input:
        train_sim_reduced_fasta=f"work_folder{pn}/subsets/train/genes/fasta/cdhit/{{selected_data}}.fasta",
        validation_sim_reduced_fasta=f"work_folder{pn}/subsets/validation/genes/fasta/cdhit/{{selected_data}}.fasta",
        test_sim_reduced_fasta=f"work_folder{pn}/subsets/test/genes/fasta/cdhit/{{selected_data}}.fasta"
    output:
        sim_train=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_train.out",
        sim_validation=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_validation.out",
        sim_test=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_test.out",
        sim_trainclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_train.out.clstr",
        sim_validationclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_validation.out.clstr",
        sim_testclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_test.out.clstr",
    log:
        f"logs{pn}/subsets/interset_similarity/cdhit/{{selected_data}}.log"
    shell:
        """
        exec > {log} 2>&1
        {params.cdhit_location}/cd-hit-2d -i {input.train_sim_reduced_fasta} -i2 {input.test_sim_reduced_fasta} \
            -o {output.sim_train} -c 0.4 -n 2 -T {threads}
        {params.cdhit_location}/cd-hit-2d -i {input.train_sim_reduced_fasta} -i2 {input.validation_sim_reduced_fasta} \
            -o {output.sim_validation} -c 0.4 -n 2 -T {threads}
        {params.cdhit_location}/cd-hit-2d -i {input.validation_sim_reduced_fasta} -i2 {input.test_sim_reduced_fasta} \
            -o {output.sim_test} -c 0.4 -n 2 -T {threads}
        """


rule get_redundant_list:
    input:
        sim_trainclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_train.out.clstr",
        sim_validationclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_validation.out.clstr",
        sim_testclr=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_test.out.clstr",
    output:
        redundant_proteins=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_redundant_proteins.txt"
    log:
        f"logs{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_redundant.log"
    shell:
        r"""
        exec > {log} 2>&1
        sed -nE 's/.*>([A-Za-z0-9-]+)....*%$/\1/p' {input.sim_trainclr} > {output.redundant_proteins}
        sed -nE 's/.*>([A-Za-z0-9-]+)....*%$/\1/p' {input.sim_validationclr} >> {output.redundant_proteins}
        sed -nE 's/.*>([A-Za-z0-9-]+)....*%$/\1/p' {input.sim_testclr} >> {output.redundant_proteins}
        """



rule cdhit_to_gene_list:
    input:
        sim_reduced_fasta=f"work_folder{pn}/subsets/{{subset}}/genes/fasta/cdhit/{{selected_data}}.fasta",
        between_reduced_redundant_proteins=f"work_folder{pn}/subsets/interset_similarity/cdhit/{{selected_data}}_redundant_proteins.txt"
    output:
        gene_list=f"work_folder{pn}/subsets/{{subset}}/genes/cdhit/genes_{{selected_data}}.txt"
    log:
        f"logs{pn}/subsets/{{subset}}/genes/cdhit/{{selected_data}}.log"
    run:
        with open(input.between_reduced_redundant_proteins, "r") as f:
            redundant_genes = set(line.strip() for line in f)

        with open(input.sim_reduced_fasta, "r") as f:
            genes = [line.strip()[1:] for line in f if line.startswith(">")]
            
        with open(output.gene_list, "w") as w:
            for gene in genes:
                if gene not in redundant_genes:
                    w.write(gene + "\n")