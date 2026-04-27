import datetime
import torch
from transformers import AutoModel, AutoTokenizer
import re
import argparse

class EmbeddWorker:
    def __init__(self, chosen_model):
        self.tokenizer, self.model = self.download_setup_model(chosen_model)

    def get_embeddings(self, sequence, i, n):
        inputs = self.tokenizer(sequence, return_tensors="pt")
        m2 = datetime.datetime.now()
        inputs = {k: v.to("cuda") for k, v in inputs.items()}
        with torch.no_grad():
            outputs = self.model(**inputs)
        emb = outputs.last_hidden_state.squeeze(0).float().cpu()
        e = datetime.datetime.now()
        print(f"Embedding_time: {e - m2} for sequence {i} / {n}")
        del inputs, outputs
        torch.cuda.empty_cache()
        return emb

    @staticmethod
    def download_setup_model(chosen_model):
        tokenizer = AutoTokenizer.from_pretrained(chosen_model)
        model = AutoModel.from_pretrained(
            chosen_model,
            torch_dtype=torch.float16,
            device_map="cuda").eval()
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


def get_all_embeddings(fasta_file, chosen_model):
    gene_name_seq_dict = read_fasta(fasta_file)
    genes = list(gene_name_seq_dict.keys())
    sequences = list(gene_name_seq_dict.values())

    em = EmbeddWorker(chosen_model)
    n = len(sequences)
    embeddings = {
        gene: em.get_embeddings(seq, i + 1, n)
        for i, (gene, seq) in enumerate(zip(genes, sequences))
    }
    return embeddings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get embeddings from protein fasta")
    parser.add_argument("--protein_fasta", required=True, help="Path to input protein fasta")
    parser.add_argument("--model_name", required=True, help="Name of embedding model (huggingface)")
    parser.add_argument("--embedding_output", required=True, help="Path to output .pt file")
    args = parser.parse_args()

    embeddings = get_all_embeddings(args.protein_fasta, args.model_name)
    torch.save(embeddings, args.embedding_output)
