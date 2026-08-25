rule get_unipro2pbd:
    #SIFTS
    output:
        unipro2pbd = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz"
    shell:
        """
        wget -O {output.unipro2pbd} https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
        """


rule download_pdb_structures:
    # flat is the union of ms+y2h, so its matched pairs cover every
    # structure any dataset could need in one shared, deduplicated cache.
    input:
        matched_pairs = "work_folder/analysis/interfaces/flat_matched_pdbs.pq"
    output:
        pdb_folder = directory("work_folder/analysis/interfaces/pdb_structures")
    log:
        "logs/analysis/interfaces/download_pdb_structures.log"
    script:
        "scripts/download_pdb_structures.py"


rule match_pairs_to_pdbs:
    input:
        pod_file = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq",
        sifts_file = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz"
    output:
        matched_pairs = "work_folder/analysis/interfaces/{dataset}_matched_pdbs.pq"
    log:
        "logs/analysis/interfaces/{dataset}_match_pairs_to_pdbs.log"
    script:
        "scripts/match_pairs_to_pdbs.py"


rule annotate_interfaces:
    input:
        matched_pairs = "work_folder/analysis/interfaces/{dataset}_matched_pdbs.pq",
        sifts_file = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz",
        pdb_folder = "work_folder/analysis/interfaces/pdb_structures"
    output:
        interface_annotated = "work_folder/analysis/interfaces/{dataset}_undirectional_interfaces.pq"
    threads: 10
    log:
        "logs/analysis/interfaces/{dataset}_annotate_interfaces.log"
    script:
        "scripts/get_interfaces.py"

rule model_interface_vs_min_size:
    input:
        gene_fasta = "work_folder/protein_sequences/uniprot_canonical.fasta",
        interface_annotated = "work_folder/analysis/interfaces/{dataset}_undirectional_interfaces.pq"
    output:
        model = "work_folder/analysis/interfaces/{dataset}_interface_size_model.joblib",
        protein_lengths = "work_folder/analysis/interfaces/{dataset}_protein_lengths.csv"
    log:
        "logs/analysis/interfaces/{dataset}_interface_size_model.log"
    script:
        "scripts/model_interface_vs_min_size.py"


rule plot_interface_size_model:
    input:
        model = "work_folder/analysis/interfaces/{dataset}_interface_size_model.joblib",
        protein_lengths = "work_folder/analysis/interfaces/{dataset}_protein_lengths.csv"
    output:
        plot = "work_folder/analysis/interfaces/{dataset}_interface_size_model.png"
    log:
        "logs/analysis/interfaces/{dataset}_plot_interface_size_model.log"
    script:
        "scripts/plot_interface_size_model.py"
