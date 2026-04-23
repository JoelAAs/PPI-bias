import torch
import pandas as pd

embeddings: dict = torch.load(snakemake.input.embeddings)

rows = []
for gene_name, emb in embeddings.items():
    mean_vec = emb.mean(dim=0).numpy()
    max_vec = emb.max(dim=0).values.numpy()
    rows.append([*mean_vec, *max_vec, gene_name])

hidden_dim = next(iter(embeddings.values())).shape[1]
columns = list(range(hidden_dim * 2)) + ["gene_name"]

pd.DataFrame(rows, columns=columns).to_csv(
    snakemake.output.protein_embeddings, sep="\t", index=False
)
