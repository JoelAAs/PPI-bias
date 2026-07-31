UNIPROT_ANNOTATION_FIELDS = (
    "accession,length,ft_transmem,cc_subcellular_location,ec,keyword,"
    "ft_signal,protein_families,cc_catalytic_activity,ft_region,lit_pubmed_id"
)
UNIPROT_COLUMN_RENAME = {
    "Entry": "uniprot_id",
    "Length": "length",
    "Transmembrane": "ft_transmem",
    "Subcellular location [CC]": "cc_subcellular_location",
    "EC number": "ec",
    "Keywords": "keyword",
    "Signal peptide": "ft_signal",
    "Protein families": "protein_families",
    "Catalytic activity": "cc_catalytic_activity",
    "Region": "ft_region",
    "PubMed ID": "lit_pubmed_id",
}





rule plot_degree_distributions:
    """
    Diagnostic for the shared protein_degrees table (module 2/3's hub bins): deg_neg
    and deg_pos distributions with the low/medium/high boundaries from
    protein_degrees.py's bin_by_quantile overlaid (taken from the bin_summary table -
    the max value observed in the "low" and "medium" bins approximates the
    lo/hi quantile cutoffs), plus a heatmap of how many proteins fall into each
    deg_neg_bin x deg_pos_bin combination (are the "high deg_neg" and "high deg_pos"
    populations mostly disjoint, as module 2's hub definition assumes?).
    """
    input:
        degree = expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree.pq",
            dataset=config["datasets"]
        ),
        bin_summary = expand(
            "work_folder/analysis/protein_degree/{dataset}_{{network_type}}_degree_bin_summary.tsv",
            dataset=config["datasets"]
        )
    output:
        distributions = "work_folder/analysis/protein_degree/plots/{network_type}_degree_distributions.png",
        heatmap = "work_folder/analysis/protein_degree/plots/{network_type}_degree_bin_heatmap.png"
    log:
        "logs/analysis/protein_degree/plots/{network_type}_degree_distributions.log"
    script:
        "scripts/plot_degree_distributions.R"
