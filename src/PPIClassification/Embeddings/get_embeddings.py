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
        token_emb = outputs.last_hidden_state.squeeze(0)  # [L, D] on GPU
        mean_vec = token_emb.mean(dim=0).float().cpu()
        max_vec = token_emb.max(dim=0).values.float().cpu()
        e = datetime.datetime.now()
        print(f"Embedding_time: {e - m2} for sequence {i} / {n}")
        del inputs, outputs, token_emb
        torch.cuda.empty_cache()
        return mean_vec, max_vec

    @staticmethod
    def download_setup_model(chosen_model):
        tokenizer = AutoTokenizer.from_pretrained(chosen_model)
        model = AutoModel.from_pretrained(chosen_model).eval()
        model = model.half().to("cuda")
        return tokenizer, model

def read_fasta(fasta_filename, accenssion=True):
    id_seq_dict = dict()
    with open(fasta_filename) as f:
        seq_id = ""
        for line in f:
            if line[0] == ">":
                if line == "":
                    continue
                if seq_id:
                    if accenssion:
                        seq_id = seq_id.groups()[0]
                    id_seq_dict[seq_id] = seq
                if accenssion:
                    seq_id = line.split("|")[1]
                else:
                    seq_id = re.search(" GN=([A-Za-z/0-9-]+) ", line)
                seq = ""
            else:
                seq += line.strip()

        if seq_id:
            seq_id = seq_id.groups()[0]
            id_seq_dict[seq_id] = seq
    return id_name_seid_seq_dictq_dict


def get_all_embeddings(fasta_file, chosen_model):
    id_name_seq_dict = read_fasta(fasta_file)
    seq_ids = list(id_seq_dict.keys())
    sequences = list(id_seq_dict.values())

    em = EmbeddWorker(chosen_model)
    n = len(sequences)
    embeddings = {
        seq_id: em.get_embeddings(seq, i + 1, n)
        for i, (seq_id, seq) in enumerate(zip(seq_ids, sequences))
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
