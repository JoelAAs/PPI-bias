import pandas as pd
import re
import networkx as nx



rule get_metis:
    params:
        script_location="src/PPIClassification/DataSplit/METIS_from_graph.py"
    input:
        graph = f"work_folder{pn}/subsets/graphs/{{graph}}.graphml"
    output:
        metis_graph = f"work_folder{pn}/subsets/graphs/metis/{{graph}}.graph",
        metis_id = f"work_folder{pn}/subsets/graphs/metis/{{graph}}_gene_id.txt"
    shell:
        """
        python3 {params.script_location} \
            --graph {input.graph} \
            --output_metis {output.metis_graph} \
            --output_int_id {output.metis_id}
        """

rule get_kahip_partitions:
    params:
        kahip_location=config["kahip_location"],
        seed=config["seed"],
        k = 20
    input:
        metis_graph=f"work_folder{pn}/subsets/graphs/metis/{{graph}}.graph"
    output:
        partitions=f"work_folder{pn}/subsets/partitions/{{graph}}.txt"
    shell:
        """
        {params.kahip_location}  {input.metis_graph} --seed={params.seed} --output_file={output.partitions} --k={params.k} --preconfiguration=strong 
        """

rule get_gene_to_partition:
    input:
        partitions = f"work_folder{pn}/subsets/partitions/{{graph}}.txt",
        gene_int_id= f"work_folder{pn}/subsets/graphs/metis/{{graph}}_gene_id.txt"
    output:
        gene_partition = f"work_folder{pn}/subsets/partitions/{{graph}}_gene_name.txt"
    run:
        rows = []
        int_id = 1
        with open(input.partitions, "r") as f:
            for line in f:
                rows.append({"int_id": int_id, "partition": line.strip()})
                int_id += 1

        partition_df = pd.DataFrame(rows)
        gene_id = pd.read_csv(input.gene_int_id, sep="\t")

        gene_partition_df = gene_id.merge(partition_df, on="int_id")
        gene_partition_df.to_csv(output.gene_partition, sep="\t", index=False)

        