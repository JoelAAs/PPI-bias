import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer
import re
import argparse

def read_fasta(fasta_filename):
    gene_name_seq_dict = dict()
    with open(fasta_filename) as f:
        gene_name = ""
        for line in f:
            if line[0] == ">":
                if gene_name:
                    gene_name_seq_dict[gene_name] = seq
                gene_name = re.search(" GN=([A-Za-z/0-9-]+) ", line).groups()[0]
                seq = ""
            else:
                seq += line.strip()
    return gene_name_seq_dict


def download_setup_model(model_name):
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, dtype=torch.float16).eval()
    return tokenizer, model


def get_mean_embeddings(fasta_file, chosen_model):
    tokenizer, model = download_setup_model(chosen_model)
    gene_name_seq_dict = read_fasta(fasta_file)
    sequences = list(gene_name_seq_dict.values())
    genes = list(gene_name_seq_dict.keys())
    inputs = tokenizer(sequences, return_tensors="pt", padding=True)
    with torch.no_grad():
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state

    mean_embeddings = embeddings.mean(dim=1)
    return mean_embeddings, genes



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--protein_fasta", required=True, help="Path to input protein fasta")
    parser.add_argument("--model_name", required=True, help="Name of embedding model (huggingface)")
    parser.add_argument("--embedding_csv", required=True, help="Path to output csv file")
    args = parser.parse_args()
    fasta_filename = args.protein_fasta
    model_name = args.model_name
    output_csv = args.embedding_csv

    mean_embeddings, genes = get_mean_embeddings(fasta_filename, model_name)
    df_embeddings = pd.DataFrame(mean_embeddings)
    df_embeddings["gene_name"] = genes
    df_embeddings.to_csv(output_csv, sep="\t", index=False)

