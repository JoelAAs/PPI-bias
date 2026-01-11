import datetime
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer
import re
import argparse

class EmbeddWorker:
    def __init__(self, chosen_model):
        self.tokenizer, self.model = self.download_setup_model(chosen_model)

    def get_mean_embeddings(self, sequences, i, n):
        inputs = self.tokenizer(sequences, return_tensors="pt", padding=True)
        m2 = datetime.datetime.now()
        inputs = {k: v.to("cuda") for k, v in inputs.items()}  # GPU transfer
        with torch.no_grad():
            outputs = self.model(**inputs)
        mean_embeddings = outputs.last_hidden_state.mean(dim=1)
        e = datetime.datetime.now()
        print(f"Embedding_time: {e - m2} for sequences {i} / {n} ")
        del inputs, outputs  # Delete tensors
        torch.cuda.empty_cache()
        return mean_embeddings

    @staticmethod
    def download_setup_model(chosen_model):
        tokenizer = AutoTokenizer.from_pretrained(chosen_model)
        model = AutoModel.from_pretrained(
            chosen_model,
            dtype=torch.float16,
            device_map="cuda",
            low_cpu_mem_usage=True).eval()
        return tokenizer, model

def read_fasta(fasta_filename):
    gene_name_seq_dict = dict()
    with open(fasta_filename) as f:
        gene_name = ""
        for line in f:
            if line[0] == ">":
                if line == "":
                    continue
                if gene_name:
                    gene_name = gene_name.groups()[0]
                    gene_name_seq_dict[gene_name] = seq
                gene_name = re.search(" GN=([A-Za-z/0-9-]+) ", line)
                seq = ""
            else:
                seq += line.strip()

        if gene_name:
            gene_name = gene_name.groups()[0]
            gene_name_seq_dict[gene_name] = seq
    return gene_name_seq_dict


def get_all_mean_embeddings(fasta_file, chosen_model, chunk_size):
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

    em = EmbeddWorker(chosen_model)
    seq_bins = binit(sequences, chunk_size)
    embeddings = [
        em.get_mean_embeddings(seqs, i*, len(sequences))
        for i, seqs in enumerate(seq_bins)
    ]
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

    chuck_size = 8
    mean_embeddings, genes = get_all_mean_embeddings(fasta_filename, model_name, chuck_size)
    df_embeddings = pd.DataFrame(mean_embeddings)
    df_embeddings["gene_name"] = genes
    df_embeddings.to_csv(output_csv, sep="\t", index=False)

