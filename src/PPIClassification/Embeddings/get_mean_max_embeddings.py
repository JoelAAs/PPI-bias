import torch
import pandas as pd

embeddings: dict = torch.load(snakemake.input.embeddings)

rows = []
for gene_name, (mean_vec, max_vec) in embeddings.items():
    rows.append([*mean_vec.numpy(), *max_vec.numpy(), gene_name])

hidden_dim = next(iter(embeddings.values()))[0].shape[0]
columns = list(range(hidden_dim * 2)) + ["gene_name"]

pd.DataFrame(rows, columns=columns).to_csv(
    snakemake.output.protein_embeddings, sep="\t", index=False
)
