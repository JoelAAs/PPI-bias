import datetime
import torch
import re
import argparse
from transformers import T5EncoderModel, T5Tokenizer
from get_embeddings import read_fasta


class EmbeddWorkerProstT5:
    def __init__(self, chosen_model):
        self.tokenizer = T5Tokenizer.from_pretrained(chosen_model, do_lower_case=False, legacy=True)
        self.model = T5EncoderModel.from_pretrained(chosen_model).eval().half().to("cuda")

    @staticmethod
    def prepare_sequence(sequence):
        # rare/ambiguous residues -> X, space-separated residues, <AA2fold> prefix for amino-acid input
        sequence = re.sub(r"[UZOB]", "X", sequence.upper())
        return "<AA2fold> " + " ".join(sequence)

    def get_embeddings(self, sequence, i, n):
        inputs = self.tokenizer(self.prepare_sequence(sequence), return_tensors="pt", add_special_tokens=True)
        m2 = datetime.datetime.now()
        inputs = {k: v.to("cuda") for k, v in inputs.items()}
        with torch.no_grad():
            outputs = self.model(**inputs)
        token_emb = outputs.last_hidden_state.squeeze(0)[1:len(sequence) + 1]  # strip <AA2fold> prefix and </s>
        mean_vec = token_emb.mean(dim=0).float().cpu()
        max_vec = token_emb.max(dim=0).values.float().cpu()
        e = datetime.datetime.now()
        print(f"Embedding_time: {e - m2} for sequence {i} / {n}")
        del inputs, outputs, token_emb
        torch.cuda.empty_cache()
        return mean_vec, max_vec


def get_all_embeddings(fasta_file, chosen_model):
    id_seq_dict = read_fasta(fasta_file)
    seq_ids = list(id_seq_dict.keys())
    sequences = list(id_seq_dict.values())

    em = EmbeddWorkerProstT5(chosen_model)
    n = len(sequences)
    embeddings = {
        seq_id: em.get_embeddings(seq, i + 1, n)
        for i, (seq_id, seq) in enumerate(zip(seq_ids, sequences))
    }
    return embeddings


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get ProstT5 embeddings from protein fasta")
    parser.add_argument("--protein_fasta", required=True, help="Path to input protein fasta")
    parser.add_argument("--model_name", required=True, help="ProstT5 model name (e.g. Rostlab/ProstT5)")
    parser.add_argument("--embedding_output", required=True, help="Path to output .pt file")
    args = parser.parse_args()

    embeddings = get_all_embeddings(args.protein_fasta, args.model_name)
    torch.save(embeddings, args.embedding_output)
