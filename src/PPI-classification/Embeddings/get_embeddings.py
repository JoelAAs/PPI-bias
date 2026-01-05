import datetime
import multiprocessing as mp
from functools import partial
import os
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
    #torch.set_num_threads(1)
    #os.environ['OMP_NUM_THREADS'] = '1'
    global tokenizer, model
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, dtype=torch.float16).eval()


def get_mean_embeddings(sequences):
    global tokenizer, model
    m1 = datetime.datetime.now()
    inputs = tokenizer(sequences, return_tensors="pt", padding=True)
    m2 = datetime.datetime.now()
    print(f"Tokenise input: {m2 - m1}")
    with torch.no_grad():
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state
    mean_embeddings = embeddings.mean(dim=1)
    e = datetime.datetime.now()
    print(f"Embedding_time: {e - m2}")

    return mean_embeddings

def get_all_mean_embeddings(fasta_file, chosen_model, chunk_size, n_cores):
    gene_name_seq_dict = read_fasta(fasta_file)
    sequences = list(gene_name_seq_dict.values())
    genes = list(gene_name_seq_dict.keys())

    def binit(x, n):
        binned = []
        while  n < len(x):
            binned.append(x[:n])
            x = x[n:]

        if len(x):
            binned.append(x)
        return binned

    seq_bins = binit(sequences, chunk_size)
    with mp.Pool(n_cores, initializer=download_setup_model, initargs=(chosen_model,)) as pool:
        embeddings = pool.map(get_mean_embeddings, seq_bins)
    embeddings = torch.cat(embeddings, dim=0)
    return embeddings, genes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get mean embeddings from protein fasta")
    parser.add_argument("--protein_fasta", required=True, help="Path to input protein fasta")
    parser.add_argument("--model_name", required=True, help="Name of embedding model (huggingface)")
    parser.add_argument("--embedding_csv", required=True, help="Path to output csv file")
    args = parser.parse_args()
    fasta_filename = args.protein_fasta
    model_name = args.model_name
    output_csv = args.embedding_csv

    chuck_size = 1000
    n_cores = 20
    mean_embeddings, genes = get_all_mean_embeddings(fasta_filename, model_name, chuck_size, n_cores)
    df_embeddings = pd.DataFrame(mean_embeddings)
    df_embeddings["gene_name"] = genes
    df_embeddings.to_csv(output_csv, sep="\t", index=False)

