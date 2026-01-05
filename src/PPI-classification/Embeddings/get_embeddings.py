import torch
from transformers import AutoModel, AutoTokenizer
import re

def read_fasta(fasta_filename):
    gene_name_seq_dict = dict()
    with open(fasta_filename) as f:
        gene_name = ""
        for line in f:
            if line[0] == ">":
                if gene_name:
                    gene_name_seq_dict[gene_name] = seq
                gene_name = re.search("GN=([A-Za-z0-9]+)\t", line).groups()[0]
                seq = ""
            else:
                seq += line.strip()
    return gene_name_seq_dict


def download_setup_model(model_name):
    #model_name = "facebook/esm2_t33_650M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, dtype=torch.float16).eval()
    return tokenizer, model


def get_mean_embeddings(fasta_file, chosen_model):
    tokenizer, model = download_setup_model(chosen_model)
    gene_name_seq_dict = read_fasta(fasta_file)
    sequences = gene_name_seq_dict.values()
    genes = gene_name_seq_dict.keys()
    inputs = tokenizer(sequences, return_tensors="pt", padding=True)
    with torch.no_grad():
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state

    mean_embeddings = embeddings.mean(dim=1)
    return mean_embeddings

model_name = "facebook/esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name, dtype=torch.float16).eval()


sequence = ["MKVLWAALLVTALAAGSLAEAAATA", "MKVLWAALLVTALAAGSLAEAAALSLSLSLSPTA"]  # Example AA sequence


pool_mean = embeddings.mean(axis=1)

