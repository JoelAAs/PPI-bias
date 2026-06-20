import pandas as pd


def main():
    similarity_df = pd.read_csv(snakemake.input.similarity_tsv, sep="\t", dtype={"qgene": "string", "sgene": "string"})
    similarity_df["edge_id"] = similarity_df[["qgene", "sgene"]].apply(
        lambda x: "-".join(sorted(x)), axis=1
    )
    similarity_df = similarity_df[["edge_id", "bitscore_p_residue"]]
    similarity_df = similarity_df.set_index("edge_id")
    with open(snakemake.output.per_edge_similarity, "w") as es_w:
        es_w.write("bait\tprey\tbitscore_p_residue\tbatch\tset\n")
        with open(snakemake.output.per_protein_similarity, "w") as ps_w:
            ps_w.write("gene\tbitscore_p_residue\tbatch\tset\n")
            for input_files, label in zip(
                [snakemake.input.perm_pos, snakemake.input.perm_neg, snakemake.input.perm_random],
                ["positive", "negative", "random_negative"]):
                for i, edge_file in enumerate(input_files):
                    edge_df = pd.read_csv(edge_file, sep="\t", dtype={"bait": "string", "prey": "string"})
                    edge_df["edge_id"] = edge_df[["bait", "prey"]].apply(
                        lambda x: "-".join(sorted(x)), axis=1
                    )
                    edge_df = edge_df.set_index("edge_id")
                    merged_df = edge_df.join(similarity_df, how="left")
                    merged_df.fillna(0, inplace=True)
                    merged_df["batch"] = i
                    merged_df["set"] = label
                    merged_df.to_csv(
                        es_w,
                        columns=["bait", "prey", "bitscore_p_residue", "batch", "set"],
                        index=False,
                        header=False,
                        sep="\t"
                    )

                    per_gene_avg = pd.concat(
                        [
                            merged_df[["bait", "bitscore_p_residue"]].rename(
                                columns={"bait": "gene"}
                            ),
                            merged_df[["prey", "bitscore_p_residue"]].rename(
                                columns={"prey": "gene"}
                            ),
                        ],
                        ignore_index=True,
                    )

                    per_gene_avg = per_gene_avg.groupby("gene").mean()
                    for gene, row in per_gene_avg.iterrows():
                        ps_w.write(f"{gene}\t{row['bitscore_p_residue']}\t{i}\t{label}\n")


if __name__ == "__main__":
    main()
