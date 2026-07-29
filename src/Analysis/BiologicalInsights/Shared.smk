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


rule download_uniprot_annotations:
    """
    UniProt REST, human reviewed proteome, bulk TSV stream. Columns renamed from
    UniProt's display headers (Entry, Length, ...) to the snake_case names used
    throughout this module. Protected download - do not re-fetch if present.
    Shared by modules 2 (hub annotation) and reference_counts below (lit_pubmed_id).
    """
    output:
        annotation=protected(config["uniprot_annotation"])
    log:
        "logs/data/uniprot/annotations.log"
    retries: 3
    run:
        import os
        import pandas as pd
        import requests

        os.makedirs(os.path.dirname(output.annotation), exist_ok=True)
        with open(log[0], "w") as logf:
            resp = requests.get(
                "https://rest.uniprot.org/uniprotkb/stream",
                params={
                    "query": "(reviewed:true) AND (organism_id:9606)",
                    "format": "tsv",
                    "fields": UNIPROT_ANNOTATION_FIELDS,
                },
            )
            if not resp.ok:
                raise RuntimeError(f"UniProt bulk download failed: HTTP {resp.status_code}")
            logf.write(f"UniProt release: {resp.headers.get('X-UniProt-Release')}, "
                       f"release date: {resp.headers.get('X-UniProt-Release-Date')}\n")
            df = pd.read_csv(pd.io.common.StringIO(resp.text), sep="\t")
            df = df.rename(columns=UNIPROT_COLUMN_RENAME)
            df.to_csv(output.annotation, sep="\t", index=False)
            logf.write(f"{len(df)} reviewed human proteins written to {output.annotation}\n")


rule reference_counts:
    """
    Per-protein publication count, used as the "study attention" proxy (module 3's
    attention confound test). Uses UniProt's own lit_pubmed_id field - the total set
    of publications UniProt has attached to that accession - rather than counting
    pubmed_ids within this PPI dataset, since the latter is circular (it only counts
    papers that already passed this project's own screening/study filters, not
    independent research attention). Shared across modules 2, 3 and 4 - built once,
    single source of truth.
    """
    input:
        uniprot_annotation = config["uniprot_annotation"]
    output:
        reference_counts = "work_folder/analysis/biological_insights/reference_counts.tsv"
    log:
        "logs/analysis/biological_insights/reference_counts.log"
    script:
        "scripts/reference_counts.py"


rule protein_degrees:
    """
    Per-protein degree table (deg_pos, deg_neg, deg_tested, sum_tests, n_pubmed,
    deg_neg_bin, deg_pos_bin). Shared by module 2 (hub definition) and module 3
    (quadrant axes) so both see identical degree numbers - built once from the POD
    table. The POD table is an already-merged undirected edge list (one row per
    unordered pair), so the two ID columns are just the two endpoints of each edge -
    which one is "bait" vs "prey" doesn't matter for degree.

    deg_pos counts any edge observed at least once, no minimum-test requirement.
    deg_neg only counts edges tested >= hub_min_neg_tested times (filtered at the
    edge level before aggregation) - a non-interaction call needs enough attempts to
    be trusted. deg_neg_bin/deg_pos_bin are the module-2 hub annotation: bin each
    column by empirical quantile of its own distribution over all proteins here (no
    normality assumption) - value above quantile(hub_quantile_cutoff) is "high",
    below quantile(1 - hub_quantile_cutoff) is "low", else "medium" - same rule both
    sides.
    """
    params:
        col_a = f"{config['id_pattern']}_bait",
        col_b = f"{config['id_pattern']}_prey",
        hub_min_neg_tested = config["hub_min_neg_tested"],
        hub_quantile_cutoff = config["hub_quantile_cutoff"]
    input:
        pod = "work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
        reference_counts = "work_folder/analysis/biological_insights/reference_counts.tsv"
    output:
        degree = "work_folder/analysis/protein_degree/{dataset}_{network_type}_degree.pq",
        bin_summary = "work_folder/analysis/protein_degree/{dataset}_{network_type}_degree_bin_summary.tsv"
    log:
        "logs/analysis/protein_degree/{dataset}_{network_type}_degree.log"
    script:
        "scripts/protein_degrees.py"


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
