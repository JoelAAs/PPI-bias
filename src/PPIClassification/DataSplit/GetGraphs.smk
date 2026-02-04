import pandas as pd
import re
import networkx as nx


def read_fasta(fasta_filename):
    gene_name_seq_dict = dict()
    with open(fasta_filename) as f:
        gene_name = ""
        for line in f:
            if line[0] == ">":
                if gene_name:
                    gene_name = gene_name.groups()[0]
                    gene_name_seq_dict[gene_name] = seq
                gene_name = re.search(" GN=([A-Za-z/0-9-]+) ",line)
                seq = ""
            else:
                seq += line.strip()

        if gene_name:
            gene_name = gene_name.groups()[0]
            gene_name_seq_dict[gene_name] = seq
    return gene_name_seq_dict


rule blast_sequence_similarity:
    params:
        n_threads=45
    input:
        fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        similarity_tsv=f"work_folder{pn}/protein_sequences/similarity/all_vs_all.tsv"
    shell:
        """
        makeblastdb -dbtype prot -in {input.fasta} -title "Gene Name SP DB"
        blastp -query {input.fasta} -db {input.fasta} \
        -outfmt "6 qseqid stitle evalue bitscore"  \
        -max_hsps 1 -num_threads {params.n_threads} -out all_vs_all.tsv
        ## Eval > 10 not reported
        """


rule get_sequence_similarity_graph:
    input:
        similarity_tsv=f"work_folder{pn}/protein_sequences/similarity/all_vs_all.tsv",
        aa_seq_fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta"
    output:
        sequence_similarity_graph=f"work_folder{pn}/subset/graphs/sequencesimilarity.graphml"
    run:
        gene_seq_dict = read_fasta(input.aa_seq_fasta)
        mean_length = round(sum([len(s) for s in gene_seq_dict.values()]) / len(gene_seq_dict))
        ava_blast_df = pd.read_csv(input.similarity_tsv,header=None,sep="\t")
        ava_blast_df.columns = ["qseqid", "stitle", "evalue", "bitscore"]
        id_dict = {
            x.split(" ")[0].split("|")[2].replace("_HUMAN",""): re.search(r' GN=([^ ]+)',x).group(1) for
            x in ava_blast_df["stitle"].unique()}

        ava_blast_df["qgene"] = ava_blast_df["qseqid"].apply(
            lambda x: id_dict[x.split("|")[2].replace("_HUMAN","")])
        ava_blast_df["sgene"] = ava_blast_df["stitle"].apply(
            lambda x: id_dict[x.split(" ")[0].split("|")[2].replace("_HUMAN","")])
        ava_blast_df["q_length"] = ava_blast_df["qgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df["s_length"] = ava_blast_df["sgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df = ava_blast_df[ava_blast_df["qgene"] != ava_blast_df["sgene"]]  # no loops
        ava_blast_df["bitscore_p_residue"] = ava_blast_df.apply(
            lambda x: round(x["bitscore"] / min([x["q_length"], x["s_length"]]) * mean_length),axis=1)
        ava_blast_df["edge_id"] = ava_blast_df[["qgene", "sgene"]].apply(lambda x: "-".join(sorted(x)),axis=1)
        ava_blast_df = ava_blast_df[~ava_blast_df["edge_id"].duplicated(keep="first")]  # NOTE: we assume B->A ~ A->B
        ava_blast_df["edge_weight"] = ava_blast_df["bitscore_p_residue"]
        G = nx.from_pandas_edgelist(ava_blast_df,"qgene","sgene",edge_attr="edge_weight")
        nx.write_graphml(G,output.sequence_similarity_graph)


rule get_min_cut_partitions:
    input:
        full_pos=f"work_folder{pn}/subsets/{{dataset}}_limit_{{pos_limit}}_pos.csv"
    output:
        ppi_graph=f"work_folder{pn}/subset/graphs/{{dataset}}_limit_{{pos_limit}}.graphml"
    run:
        pos_df = pd.read_csv(input.full_pos,sep="\t")
        pos_df["edge_weight"] = 1
        G = nx.from_pandas_edgelist(ava_blast_df,"gene_name_bait","gene_name_prey",edge_attr="edge_weight")
        nx.write_graphml(G,output.ppi_graph)
