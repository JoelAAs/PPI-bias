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


if __name__ == "__main__":
    localisation_file = snakemake.input.annotation
    gene_names_file = snakemake.input.gene_names
    id_pattern = snakemake.params.id_pattern
    output_file = snakemake.output.protein_localisation

    df_localisation = pd.read_csv(localisation_file, sep="\t")
    df_localisation = df_localisation[
        df_localisation["Reliability"].isin(["Supported", "Approved", "Enhanced"])
    ]

    ensembl_entrez = ensembl_to_entrez(df_localisation["Gene"].tolist())

    rows = []
    for _, row in df_localisation.iterrows():
        entrez_id = ensembl_entrez.get(row["Gene"])
        if entrez_id is None:
            continue
        for loc in str(row["Main location"]).split(";"):
            loc = loc.strip()
            if loc and loc != "nan":
                rows.append((entrez_id, loc))
    df_entrez_localisation = pd.DataFrame(rows, columns=["gene_name", "localisation"]).drop_duplicates()

    if id_pattern == "uniprot_id":
        df_gene_names = pd.read_csv(gene_names_file, sep="\t", dtype=str)
        df_out = df_gene_names.merge(df_entrez_localisation, on="gene_name")[["uniprot_id", "localisation"]]
        df_out = df_out.rename(columns={"uniprot_id": "id"})
    else:
        df_out = df_entrez_localisation.rename(columns={"gene_name": "id"})

    df_out.drop_duplicates().to_csv(output_file, sep="\t", index=False)
