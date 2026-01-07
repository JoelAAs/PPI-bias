import datetime
import ray
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer
import re
import argparse

@ray.remote(num_cpus=1)
class RayEmbeddWorker:
    def __init__(self, chosen_model):
        self.tokenizer, self.model = self.download_setup_model(chosen_model)

    def get_mean_embeddings(self, sequences, i, n):
        inputs = self.tokenizer(sequences, return_tensors="pt", padding=True)
        m2 = datetime.datetime.now()
        with torch.no_grad():
            outputs = self.model(**inputs)
            embeddings = outputs.last_hidden_state
        mean_embeddings = embeddings.mean(dim=1)
        e = datetime.datetime.now()
        print(f"Embedding_time: {e - m2} for sequences {i-len(sequences)} - {i} / {n} ")

        return mean_embeddings

    @staticmethod
    def download_setup_model(chosen_model):
        tokenizer = AutoTokenizer.from_pretrained(chosen_model)
        model = AutoModel.from_pretrained(chosen_model, dtype=torch.float16).eval()
        return tokenizer, model

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


    ray.init(num_cpus=n_cores)

    workers = [RayEmbeddWorker.remote(chosen_model) for _ in range(n_cores)]
    print("Ray workers recruited")
    seq_bins = binit(sequences, chunk_size)
    work_queue = []
    for i, seqs in enumerate(seq_bins):
        chosen_worker = workers[i%len(workers)]
        work_queue.append(chosen_worker.get_mean_embeddings.remote(seqs, i*chunk_size, len(sequences)))

    embeddings = ray.get(work_queue)
    embeddings = torch.cat(embeddings, dim=0)
    ray.shutdown()

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

    chuck_size = 5
    n_cores = 40
    mean_embeddings, genes = get_all_mean_embeddings(fasta_filename, model_name, chuck_size, n_cores)
    df_embeddings = pd.DataFrame(mean_embeddings)
    df_embeddings["gene_name"] = genes
    df_embeddings.to_csv(output_csv, sep="\t", index=False)

