import datetime
import torch
import re
import argparse
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig


class EmbeddWorkerESMC:
    def __init__(self, model_name):
        self.model = ESMC.from_pretrained(model_name).to("cuda").eval()

    def get_embeddings(self, sequence, i, n):
        t0 = datetime.datetime.now()
        protein = ESMProtein(sequence=sequence)
        protein_tensor = self.model.encode(protein)
        with torch.no_grad():
            output = self.model.logits(protein_tensor, LogitsConfig(embeddings=True))
        emb = output.embeddings[1:-1].float().cpu()  # strip BOS/EOS tokens
        print(f"Embedding_time: {datetime.datetime.now() - t0} for sequence {i} / {n}")
        torch.cuda.empty_cache()
        return emb


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


def get_all_embeddings(fasta_file, model_name):
    gene_name_seq_dict = read_fasta(fasta_file)
    genes = list(gene_name_seq_dict.keys())
    sequences = list(gene_name_seq_dict.values())

    em = EmbeddWorkerESMC(model_name)
    n = len(sequences)
    return {
        gene: em.get_embeddings(seq, i + 1, n)
        for i, (gene, seq) in enumerate(zip(genes, sequences))
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get ESMC embeddings from protein fasta")
    parser.add_argument("--protein_fasta", required=True, help="Path to input protein fasta")
    parser.add_argument("--model_name", required=True, help="ESMC model name (e.g. esmc_600m)")
    parser.add_argument("--embedding_output", required=True, help="Path to output .pt file")
    args = parser.parse_args()

    embeddings = get_all_embeddings(args.protein_fasta, args.model_name)
    torch.save(embeddings, args.embedding_output)
