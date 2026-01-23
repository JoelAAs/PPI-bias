import pandas as pd
import re
import networkx as nx


#TODO DRY fix
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


rule get_METIS_adjacency_list:
    params:
        normalize=True
    input:
        similarity_tsv=f"work_folder{pn}/protein_sequences/similarity/all_vs_all.tsv",
        aa_seq_fasta=f"work_folder{pn}/protein_sequences/gene_name_sp_dedub.fasta",
        full_pos= f"work_folder{pn}/subsets/{{dataset}}_full_limit_poslimit_{{pos_limit}}_pos.csv"
    output:
        similarity_edge_list=f"work_folder{pn}/protein_sequences/similarity/{{dataset}}/avg_bitscore_all_poslimit_{{pos_limit}}.tsv",
        gene_int_id = f"work_folder{pn}/protein_sequences/similarity/{{dataset}}/gene_int_id_poslimit_{{pos_limit}}.tsv",
        similarity_mentis=f"work_folder{pn}/protein_sequences/similarity/{{dataset}}/avg_bitscore_poslimit_{{pos_limit}}.graph"
    run:
        # Node weights
        pos_edges = pd.read_csv(input.full_pos, sep="\t")
        n_pos_edges = pos_edges.shape[0]
        all_genes_in_pos = set(pos_edges["gene_name_bait"]) | set(pos_edges["gene_name_prey"])
        bait_count =  pos_edges.groupby("gene_name_bait", as_index=False).size().rename(
            {"gene_name_bait":"gene_name"}, axis=1)
        prey_count = pos_edges.groupby("gene_name_prey",as_index=False).size().rename(
            {"gene_name_prey": "gene_name"},axis=1)
        edge_per_gene = bait_count.merge(prey_count, on="gene_name", how="outer").fillna(0) # NOTE: assumes no homodimers
        edge_per_gene["frequency"] = (edge_per_gene["size_x"] + edge_per_gene["size_y"])/n_pos_edges
        node_weights = {g:w for i, (g, w) in edge_per_gene[["gene_name", "frequency"]].iterrows()}

        #
        gene_seq_dict = read_fasta(input.aa_seq_fasta)
        mean_length = round(sum([len(s) for s in gene_seq_dict.values()]) / len(gene_seq_dict))
        ava_blast_df = pd.read_csv(input.similarity_tsv,header=None,sep="\t")
        ava_blast_df.columns = ["qseqid", "stitle", "evalue", "bitscore"]
        id_dict = {
            x.split(" ")[0].split("|")[2].replace("_HUMAN",""): re.search(r' GN=([^ ]+)',x).group(1) for
            x in ava_blast_df["stitle"].unique()
        }

        ava_blast_df["qgene"] = ava_blast_df["qseqid"].apply(lambda x: id_dict[x.split("|")[2].replace("_HUMAN","")])
        ava_blast_df["sgene"] = ava_blast_df["stitle"].apply(lambda x: id_dict[
            x.split(" ")[0].split("|")[2].replace("_HUMAN","")])
        ava_blast_df["q_length"] = ava_blast_df["qgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df["s_length"] = ava_blast_df["sgene"].apply(lambda x: len(gene_seq_dict.get(x,"")))
        ava_blast_df = ava_blast_df[ava_blast_df["qgene"] != ava_blast_df["sgene"]]  # no loops
        ava_blast_df["bitscore_p_residue"] = ava_blast_df.apply(
            lambda x: round(x["bitscore"] / min([x["q_length"], x["s_length"]]) * mean_length),axis=1)
        ava_blast_df["edge_id"] = ava_blast_df[["qgene", "sgene"]].apply(lambda x: "-".join(sorted(x)),axis=1)
        ava_blast_df = ava_blast_df[~ava_blast_df["edge_id"].duplicated(keep="first")] # NOTE: we assume B->A ~ A->B

        G = nx.from_pandas_edgelist(ava_blast_df,"qgene","sgene",edge_attr="bitscore_p_residue")
        sorted_nodes = list(enumerate(sorted(G.nodes())))
        int_mapping = {node: i+1 for i, node in sorted_nodes} # +1 for 1 index
        with open(output.gene_int_id, "w") as w:
            w.write("gene_name\tint_id\n")
            for gene, int_id in int_mapping.items():
                w.write(f"{gene}\t{int_id}\n")

        ava_blast_df[["qgene", "sgene", "bitscore_p_residue"]].to_csv(output.similarity_edge_list)

        with open(output.similarity_mentis,"w") as w:
            w.write(f'{G.number_of_nodes()} {G.number_of_edges()} 11\n') # 11 for node and edge weights
            for i, node in sorted_nodes:
                line = f"{node_weights[node]}  "
                line += " ".join([
                    f'{int_mapping[edge[1]]} {edge[2]["bitscore_p_residue"]}' for
                    edge in G.edges(node,data=True)])
                w.write(line + "\n")


rule get_kahip_partitions:
    params:
        kahip_location=config["kahip_location"],
        seed=config["seed"],
        k = 12
    input:
        similarity_mentis=f"work_folder{pn}/protein_sequences/similarity/avg_bitscore_poslimit_{{pos_limit}}.graph"
    output:
        partitions=f"work_folder{pn}/protein_sequences/similarity/partitions/clusters_kaffpa_poslimit_{{pos_limit}}.txt"
    shell:
        """
        {params.kahip_location}  {input.similarity_mentis} --seed={params.seed} --output_file={output.partitions} --k={params.k} --preconfiguration=strong 
        """

rule get_gene_to_partition:
    input:
        partitions = f"work_folder{pn}/protein_sequences/similarity/partitons/clusters_kaffpa_poslimit_{{pos_limit}}.txt",
        gene_int_id= f"work_folder{pn}/protein_sequences/similarity/gene_int_id_poslimit_{{pos_limit}}.tsv"
    output:
        gene_partition = f"work_folder{pn}/protein_sequences/similarity/gene_partition_poslimit_{{pos_limit}}.tsv"
    run:
        rows = []
        int_id = 1
        with open(input.partitions, "r") as f:
            for line in f:
                rows.append({"int_id": int_id, "sequence_partition": line.strip()})
                int_id += 1

        partition_df = pd.DataFrame(rows)
        gene_id = pd.read_csv(input.gene_int_id, sep="\t")

        gene_partition_df = gene_id.merge(partition_df, on="int_id")
        gene_partition_df.to_csv(output.gene_partition, sep="\t", index=False)