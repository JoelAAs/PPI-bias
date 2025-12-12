from collections import defaultdict
from mygene import MyGeneInfo
import pandas as pd


def get_go_genes(genes):
    mg = MyGeneInfo()
    result = mg.querymany(
        genes,
        scopes="symbol,alias",
        fields="go",
        species="human",
        returnall=True
    )

    gene_go_dict = dict()
    for go_q in result["out"]:
        gene_gos = dict()
        gene = go_q["query"]
        go_terms = go_q.get('go',{})
        for i, term in enumerate(["BP", "MF", "CC"]):
            go_match = go_terms.get(term,{})
            if isinstance(go_match,list):
                gos = {go["id"] for go in go_match}
            elif "id" in go_match:
                gos = {go_match["id"], }
            else:
                gos = set()
            gene_gos[term] = gos
        gene_go_dict[gene] = gene_gos

    for missing_gene in result["missing"]:
        gene_gos = {
            "BP": set(),
            "CC": set(),
            "MF": set()
        }
        gene_go_dict[missing_gene] = gene_gos
    return gene_go_dict

def get_go_frequency(go_dict, category="BP"):
    n_genes = len(go_dict)
    go_freq_dict = defaultdict(int)
    for gene in go_dict:
        for go_term in go_dict[gene][category]:
            go_freq_dict[go_term] += 1

    go_freq_df = pd.DataFrame(go_freq_dict.items(), columns = ["go_term", "go_occurrence"])
    go_freq_df["go_frequency"] = go_freq_df["go_occurrence"]/n_genes
    return go_freq_df



