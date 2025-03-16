import pandas as pd

rule placeholder_POD_calc:
    input:
        uniprot_to_gene = "work_folder/intact/uniprot_to_gene_name.csv"
    output:
        prey_pod = "work_folder/cell_type_pod/{cell_line}_pod.csv",
    run:
        pod_file = config["cell_lines"][wildcards.cell_line]["pod"]

        if pod_file:
            pod_df = pd.read_csv(pod_file, sep="\t")
            uniprot_gene_name_df = pd.read_csv(input.uniprot_to_gene, sep = "\t")
            pod_df = pod_df.merge(uniprot_gene_name_df, on="gene_name", how="inner")
            pod_df.to_csv(
                output.prey_pod,
                sep="\t",
                index=False,
                columns = [
                    "gene_name",
                    "uniprot_id",
                    "relative_frequence"
                ]

            )
        else:
            with open(output.prey_pod, "w") as w:
                w.write("gene_name\tuniprot_id\trelative_frequence\n")
            shell(f"touch {output.prey_pod}")