import pandas as pd
from mygene import MyGeneInfo


def ensembl_to_entrez(ensembl_ids):
    mg = MyGeneInfo()
    result = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="entrezgene",
        species="human",
        returnall=True
    )
    mapping = {}
    for r in result["out"]:
        entrez = r.get("entrezgene")
        if entrez is not None:
            mapping[r["query"]] = str(int(entrez))
    return mapping


def get_go_genes(entrez_ids):
    mg = MyGeneInfo()
    result = mg.querymany(
        entrez_ids,
        scopes="entrezgene",
        fields="go",
        species="human",
        returnall=True
    )

    rows = []
    for go_q in result["out"]:
        gene = go_q["query"]
        go_terms = go_q.get('go', {})
        go_match = go_terms.get("BP", {})
        if isinstance(go_match, list):
            gos = {go["id"] for go in go_match}
        elif "id" in go_match:
            gos = {go_match["id"]}
        else:
            gos = set()
        for go in gos:
            rows.append((gene, go))

    return pd.DataFrame(rows, columns=["gene_id", "go_id"])


if __name__ == "__main__":
    gene_names_file = snakemake.input.gene_names
    localisation_input = snakemake.input.annotation
    annotation_file = snakemake.output.annotation_proteins

    # gene_names_file: uniprot_id \t entrez_id (header row)
    df_gene_names = pd.read_csv(gene_names_file, sep="\t")
    entrez_ids = set(df_gene_names["gene_name"].dropna().astype(str).tolist())

    # Read localisation; use only main localisation, supported, approved or enhanced
    df_localisation = pd.read_csv(localisation_input, sep="\t")
    df_localisation = df_localisation[
        df_localisation["Reliability"].isin(["Supported", "Approved", "Enhanced"])
    ]

    # Convert Ensembl IDs (Gene column) to Entrez IDs
    ensembl_ids = df_localisation["Gene"].tolist()
    ensembl_entrez = ensembl_to_entrez(ensembl_ids)

    localisation_rows = []
    for _, row in df_localisation.iterrows():
        entrez_id = ensembl_entrez.get(row["Gene"])
        if entrez_id is None:
            continue
        for loc in str(row["Main location"]).split(";"):
            loc = loc.strip()
            if loc and loc != "nan":
                localisation_rows.append((entrez_id, loc))
    df_loc_pairs = pd.DataFrame(localisation_rows, columns=["gene_id", "annotation"])

    # Get GO biological process annotations for all genes in dataset
    df_go = get_go_genes(entrez_ids).rename(columns={"go_id": "annotation"})

    # Write localisation and go to the same file; keep only annotations with >400 genes
    df_combined = pd.concat([df_loc_pairs, df_go], ignore_index=True)
    annotation_gene_counts = df_combined.groupby("annotation")["gene_id"].nunique()
    keep = annotation_gene_counts[annotation_gene_counts > 400].index
    df_combined[df_combined["annotation"].isin(keep)].to_csv(
        annotation_file, sep="\t", index=False, header=False
    )
