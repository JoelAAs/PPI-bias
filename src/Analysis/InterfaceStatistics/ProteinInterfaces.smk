rule get_unipro2pbd:
    #SIFTS
    output:
        unipro2pbd = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz"
    shell:
        """
        wget -O {output.unipro2pbd} https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
        """


rule download_pdb_structures:
    input:
        matched_pairs = "work_folder/analysis/interfaces/flat_matched_pdbs.pq"
    output:
        pdb_folder = directory("work_folder/analysis/interfaces/pdb_structures")
    log:
        "logs/analysis/interfaces/download_pdb_structures.log"
    script:
        "scripts/download_pdb_structures.py"


rule match_pairs_to_pdbs:
    "Flat is super set of y2h/ms"
    input:
        pod_file = "work_folder/analysis/POD/undirectional/POD_flat.pq",
        sifts_file = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz"
    output:
        matched_pairs = "work_folder/analysis/interfaces/all_matched_pdbs.pq"
    log:
        "logs/analysis/interfaces/match_pairs_to_pdbs.log"
    script:
        "scripts/match_pairs_to_pdbs.py"


rule subset_matched_pdbs:
    input:
        all_matched_pairs = "work_folder/analysis/interfaces/all_matched_pdbs.pq",
        pod_file = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq"
    output:
        matched_pairs = "work_folder/analysis/interfaces/{dataset}_matched_pdbs.pq"
    log:
        "logs/analysis/interfaces/{dataset}_subset_matched_pdbs.log"
    script:
        "scripts/subset_matched_pdbs.py"


rule set_interface_size:
    input:
        matched_pairs = "work_folder/analysis/interfaces/{dataset}_matched_pdbs.pq",
        sifts_file = "work_folder/analysis/interface/pdb_chain_uniprot.tsv.gz",
        pdb_folder = "work_folder/analysis/interfaces/pdb_structures",
        gene_fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
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
        model = "work_folder/analysis/interfaces/model/{dataset}_interface_size_model.joblib",
        protein_lengths = "work_folder/analysis/interfaces/{dataset}_protein_lengths.csv"
    log:
        "logs/analysis/interfaces/{dataset}_interface_size_model.log"
    script:
        "scripts/model_interface_vs_min_size.py"


rule plot_interface_size_model:
    input:
        model = "work_folder/analysis/interfaces/model/{dataset}_interface_size_model.joblib",
        protein_lengths = "work_folder/analysis/interfaces/{dataset}_protein_lengths.csv"
    output:
        plot = "work_folder/analysis/interfaces/plots/{dataset}_interface_size_model.png"
    log:
        "logs/analysis/interfaces/{dataset}_plot_interface_size_model.log"
    script:
        "scripts/plot_interface_size_model.py"


rule plot_detection_vs_interface_size:
    input:
        interface_annotated = "work_folder/analysis/interfaces/{dataset}_undirectional_interfaces.pq"
    output:
        plot = "work_folder/analysis/interfaces/plots/{dataset}_detection_vs_interface_size.png"
    log:
        "logs/analysis/interfaces/{dataset}_plot_detection_vs_interface_size.log"
    script:
        "scripts/plot_detection_vs_interface_size.py"



rule model_detection_vs_interface:
    input:
        interface_annotated = "work_folder/analysis/interfaces/{dataset}_undirectional_interfaces.pq",
        gene_fasta = "work_folder/protein_sequences/uniprot_canonical.fasta"
    params:
        n_bootstrap = 500,
        n_permutations = 200
    output:
        coefficients = "work_folder/analysis/interfaces/{dataset}_detection_interface_coefs.tsv",
        resamples = "work_folder/analysis/interfaces/{dataset}_detection_interface_resamples.npz",
        partial_effect = "work_folder/analysis/interfaces/plots/{dataset}_detection_interface_effect.png"
    log:
        "logs/analysis/interfaces/{dataset}_detection_vs_interface.log"
    script:
        "scripts/model_detection_ratio.py"


rule plot_or_interface:
    input:
        coefficients = "work_folder/analysis/interfaces/{dataset}_detection_interface_coefs.tsv",
    output:
        or_png = "work_folder/analysis/interfaces/plots/{dataset}_detection_interface_OR.png"
    run:

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        terms = ["log_interface", "log_min_length", "log_n_tested"]
        labels = ["Interface size", "Minimum protein length", "Number tested"]

        coefs = pd.read_csv(input.coefficients, sep = "\t")
        d = coefs[coefs.term.isin(terms)].copy()
        d["label"] = d["term"].map(dict(zip(terms, labels)))
        d = d.set_index("term").loc[terms].reset_index()

        y = np.arange(len(d))
        or_ = d["or_per_doubling"]
        lo = d["or_bootstrap_low"]
        hi = d["or_bootstrap_high"]

        fig, ax = plt.subplots(figsize=(6, 3.5))

        ax.errorbar(
            or_, y,
            xerr=[or_ - lo, hi - or_],
            fmt="o", capsize=4
        )

        ax.axvline(1, ls="--", color="grey")
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels(d["label"])
        ax.invert_yaxis()
        ax.set_xlabel("Odds ratio per doubling")

        for i, p in enumerate(d["p_value_naive"]):
            ax.text(1.02, i, f"p = {p:.3g}",
                    transform=ax.get_yaxis_transform(), va="center")

        fig.tight_layout()
        fig.savefig(output.or_png, dpi=300, bbox_inches="tight")
        plt.close(fig)