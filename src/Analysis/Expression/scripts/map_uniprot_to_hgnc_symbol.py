import pandas as pd
from mygene import MyGeneInfo

if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    df = pd.read_csv(snakemake.input.gene_names, sep="\t")
    n_in = len(df)

    df["gene_name"] = df["gene_name"].astype(str)
    entrez_ids = df["gene_name"].dropna().unique().tolist()

    mg = MyGeneInfo()
    result = mg.querymany(
        entrez_ids, scopes="entrezgene", fields="symbol", species="human", returnall=True
    )
    entrez_to_symbol = {r["query"]: r["symbol"] for r in result["out"] if "symbol" in r}

    df["hgnc_symbol"] = df["gene_name"].map(entrez_to_symbol)
    n_mapped = df["hgnc_symbol"].notna().sum()
    print(f"{n_in} uniprot IDs in, {len(entrez_ids)} distinct Entrez IDs queried, "
          f"{n_mapped} uniprot IDs mapped to an HGNC symbol, "
          f"{n_in - n_mapped} dropped (unmapped, not imputed)", file=log, flush=True)

    df[["uniprot_id", "hgnc_symbol"]].dropna().to_csv(
        snakemake.output.symbol_map, sep="\t", index=False
    )
    log.close()
