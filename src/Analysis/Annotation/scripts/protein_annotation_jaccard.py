import pandas as pd
from mygene import MyGeneInfo


def map_uniprot_to_gene_id(uniprot_file):
    # gene_names on format: uniprot_id \t gene_id, with header
    map_dict = {}
    with open(uniprot_file, 'r') as f:
        next(f)
        for line in f:
            uniprot_id, gene_id = line.strip().split('\t')
            map_dict[uniprot_id] = gene_id
    return map_dict


def get_go_bp_terms(entrez_ids):
    mg = MyGeneInfo()
    result = mg.querymany(
        entrez_ids,
        scopes="entrezgene",
        fields="go",
        species="human",
        returnall=True
    )

    go_terms = {}
    for go_q in result["out"]:
        gene = go_q["query"]
        go_bp = go_q.get("go", {}).get("BP", {})
        if isinstance(go_bp, list):
            gos = {go["id"] for go in go_bp}
        elif "id" in go_bp:
            gos = {go_bp["id"]}
        else:
            gos = set()
        go_terms[gene] = gos
    return go_terms


def jaccard(set_a, set_b):
    union = set_a | set_b
    if not union:
        return float("nan")
    return len(set_a & set_b) / len(union)


if __name__ == "__main__":
    pod_file = snakemake.input.pod
    gene_names_file = snakemake.input.gene_names
    jaccard_file = snakemake.output.jaccard

    bait_column = f"{snakemake.params.id_pattern}_bait"
    prey_column = f"{snakemake.params.id_pattern}_prey"

    df_pod = pd.read_parquet(pod_file)
    df_pod = df_pod[df_pod["n_observed"] > 0].copy()

    proteins = pd.unique(df_pod[[bait_column, prey_column]].to_numpy().ravel())

    uniprot_to_gene_id = map_uniprot_to_gene_id(gene_names_file)
    protein_to_entrez = {p: uniprot_to_gene_id[p] for p in proteins if p in uniprot_to_gene_id}

    go_by_entrez = get_go_bp_terms(set(protein_to_entrez.values()))
    protein_go_terms = {
        p: go_by_entrez.get(entrez_id, set())
        for p, entrez_id in protein_to_entrez.items()
    }

    df_pod["jaccard"] = [
        jaccard(protein_go_terms.get(bait, set()), protein_go_terms.get(prey, set()))
        for bait, prey in zip(df_pod[bait_column], df_pod[prey_column])
    ]

    df_pod[[bait_column, prey_column, "lower_bound_pod", "jaccard"]].to_csv(jaccard_file, sep="\t", index=False)
