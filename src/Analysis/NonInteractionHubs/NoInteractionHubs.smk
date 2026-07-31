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
    output:
        annotation="work_folder/data/uniprot/human_annotations.tsv"
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
                    "query": "organism_id:9606",
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


rule reference_counts:
    input:
        uniprot_annotation = "work_folder/data/uniprot/human_annotations.tsv"
    output:
        reference_counts = "work_folder/analysis/biological_insights/reference_counts.tsv"
    log:
        "logs/analysis/biological_insights/reference_counts.log"
    script:
        "scripts/reference_counts.py"


rule protein_degrees:
    params:
        col_a = f"{config['id_pattern']}_bait",
        col_b = f"{config['id_pattern']}_prey",
        interaction_hubs=config["hub_limits"],
        non_interaction_hub_q = config["non_interaction_hub_quantile_cutoff"]
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


rule define_hub_set:
    input:
        degree="work_folder/analysis/protein_degree/{dataset}_{network_type}_degree.pq"
    output:
        int_hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.log"
    run:
        df = pd.read_parquet(input.degree)
        df["hub_status"] = np.where(
            (df["deg_neg_bin"] == "low") & (df["deg_pos_bin"] == "high"), "interaction",
            np.where((df["deg_neg_bin"] == "high") & (df["deg_pos_bin"] == "low"), "non-interaction", "")
        )
        df[df["hub_status"] != ""].to_csv(output.int_hub_set, sep="\t", index=None)

        
rule fetch_disorder_content:
    input:
        hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv",
        uniprot_annotation="work_folder/data/uniprot/human_annotations.tsv"
    output:
        disorder="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.log"
    script:
        "scripts/fetch_disorder_content.py"


rule test_hub_properties:
    params:
        bootstrap_replicates=config["bootstrap_replicates"],
        seed=config["seed"]
    input:
        hub_set="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_hub_set.tsv",
        uniprot_annotation="work_folder/data/uniprot/human_annotations.tsv",
        disorder="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_disorder.tsv"
    output:
        enrichment="work_folder/analysis/no_interaction_hubs/{dataset}_{network_type}_enrichment.tsv"
    log:
        "logs/analysis/no_interaction_hubs/{dataset}_{network_type}_enrichment.log"
    script:
        "scripts/test_hub_properties.py"


rule plot_hub_properties:
    input:
        expand(
            "work_folder/analysis/no_interaction_hubs/{dataset}_{{network_type}}_enrichment.tsv",
            dataset=config["datasets"]
        )
    output:
        binary="work_folder/analysis/no_interaction_hubs/plots/{network_type}_hub_properties_binary.png",
        continuous="work_folder/analysis/no_interaction_hubs/plots/{network_type}_hub_properties_continuous.png"
    log:
        "logs/analysis/no_interaction_hubs/plots/{network_type}_hub_properties.log"
    script:
        "scripts/plots/plot_hub_properties.R"
