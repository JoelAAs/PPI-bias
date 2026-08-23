rule get_unipro2pbd:
    #SIFTS
    output:
        unipro2pbd = "results/analysis/interface/pdb_chain_uniprot.tsv.gz"
    shell:
        """
        wget -O {output.unipro2pbd} https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
        """


rule annotate_interfaces:
    input:
        pod_file = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq",
        sifts_file = "results/analysis/interface/pdb_chain_uniprot.tsv.gz"
    output:
        pdb_folder = directory("work_folder/analysis/interfaces/{dataset}_pdb_structures"),
        interface_annotated = "work_folder/analysis/interfaces/{dataset}_undirectional_interfaces.pq"
    log:
        "logs/analysis/interfaces/{dataset}_annotate_interfaces.log"
    script:
        "scripts/get_interfaces.py"