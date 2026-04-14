import pandas as pd
import re


def parse_cellosaurus(filename, output_file):
    header_passed = False
    headers_to_keep = ["AC", "OX"]
    with open(output_file, "w") as w:
        w.write("\t".join(headers_to_keep) +"\n")
        with open(filename, "r") as f:
            entry = {h:[] for h in headers_to_keep}
            for line in f:
                if header_passed:
                    if line[:2] == "//":
                        line_out = "\t".join([
                            "<->".join(entry[head]) for head in entry
                        ]) + "\n"
                        w.write(line_out)
                        entry = {h:[] for h in headers_to_keep}
                        continue
                    header, value = line.strip().split("   ")
                    if header in entry:
                        entry[header].append(value)
                if line == "____________________________________________________________________________\n":
                    header_passed = True


def format_bioplex(bp_f, gene_name_df, cl):
        bp_df = pd.read_csv(bp_f, sep="\t")[["Bait Symbol", "Prey Symbol"]]
        bp_df.columns = ["gene_name_bait", "gene_name_prey"]
        bp_df[["pubmed_id", "detection_method", "CVCL"]] = ["33961781", "MS-0006", cl]


        bp_df = bp_df.merge(gene_name_df, left_on="gene_name_bait", right_on="gene_name")
        del bp_df["gene_name"]
        bp_df = bp_df.merge(gene_name_df, left_on="gene_name_prey", right_on="gene_name", suffixes=("_bait", "_prey"))
        del bp_df["gene_name"]

        return bp_df



rule get_pubtator:
    output:
        f"work_folder{pn}/cell_line_annotation/cellline2pubtator3.txt"
    log:
        f"logs{pn}/cell_line_annotation/cellline2pubtator3.log"
    shell:
        """
        exec > {log} 2>&1
        wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/cellline2pubtator3.gz
        gunzip cellline2pubtator3.gz
        mv cellline2pubtator3 {output}
        """


rule get_bioplex:
    output:
        cvcl_0291_bp=f"work_folder{pn}/data/bioplex/CVCL_0291.csv",
        cvcl_0063_bp=f"work_folder{pn}/data/bioplex/CVCL_0063.csv"
    log:
        f"logs{pn}/data/bioplex/bioplex.log"
    shell:
        """
        exec > {log} 2>&1
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_293T_DirectedEdges.tsv -O {output.cvcl_0063_bp}
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_HCT116_DirectedEdges.tsv -O {output.cvcl_0291_bp}
        """

rule get_cellosaurus:
    output:
        cellosaurus = f"work_folder{pn}/data/cellosaurus/cellosaurus.txt"
    log:
        f"logs{pn}/data/cellosaurus/cellosaurus.log"
    shell:
        """
        wget https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt -O {output.cellosaurus} > {log} 2>&1
        """

rule get_cellosaurus_human_cl:
    input:
        cellosaurus = f"work_folder{pn}/data/cellosaurus/cellosaurus.txt"
    output:
        cellosaurus_csv = f"work_folder{pn}/data/cellosaurus/taxon_cellosaurus.csv",
        cellosaurus_human = f"work_folder{pn}/data/cellosaurus/cellosaurus_human.txt"
    log:
        f"logs{pn}/data/cellosaurus/cellosaurus_human.log"
    run:
        parse_cellosaurus(input.cellosaurus, output.cellosaurus_csv)
        df_cellosaurus = pd.read_csv(output.cellosaurus_csv, sep="\t")
        df_cellosaurus = df_cellosaurus[df_cellosaurus["OX"].str.contains("9606")]
        df_cellosaurus["AC"].to_csv(output.cellosaurus_human, index=None)




rule join_to_pid:
    params:
        pid_to_drop = config["pid_to_remove"]
    input:
        pubator = f"work_folder{pn}/cell_line_annotation/cellline2pubtator3.txt",
        human_cell_lines = f"work_folder{pn}/data/cellosaurus/cellosaurus_human.txt",
        formatted_intact = f"work_folder{pn}/formated/bait_prey_publications.csv",
        cvcl_0291_bp=f"work_folder{pn}/data/bioplex/CVCL_0291.csv",
        cvcl_0063_bp=f"work_folder{pn}/data/bioplex/CVCL_0063.csv",
        manual_curated = "data/checked_studies.csv",
        gene_names = f"work_folder{pn}/gene_names/gene_names.csv"
    output:
        cvcl_ppi = f"work_folder{pn}/formated/bait_prey_CVCL.csv"
    log:
        f"logs{pn}/formated/bait_prey_CVCL.log"
    run:
        column_order =[
            'uniprot_id_bait', 'uniprot_id_prey', 'pubmed_id', 'detection_method',
            'gene_name_bait', 'gene_name_prey', 'CVCL'
        ]
        human_cl = pd.read_csv(input.human_cell_lines)

        # CL annotation
        pubtator3 = pd.read_csv(
            input.pubator,
            sep="\t",
            header=None)
        pubtator3.columns = ["pubmed_id", "annotation", "CVCL", "cell_name", "source"]
        pubtator3 = pubtator3[pubtator3["CVCL"].isin(human_cl["AC"])]
        n_cl_pid = pubtator3.groupby("pubmed_id", as_index=False)["CVCL"].nunique()
        single_pid = n_cl_pid[n_cl_pid["CVCL"] == 1]["pubmed_id"]
        single_pid_df = pubtator3[pubtator3["pubmed_id"].isin(single_pid)][["pubmed_id", "CVCL"]]
        single_pid_df["pubmed_id"] = single_pid_df.pubmed_id.astype(str)
        manual_df = pd.read_csv(input.manual_curated,sep="\t",dtype={"pubmed_id": object})[["pubmed_id","cl_id"]]
        manual_df.columns = ["pubmed_id", "CVCL"]
        manual_df = manual_df.dropna()
        
        intact_df = pd.read_csv(input.formatted_intact,sep="\t",dtype={"pubmed_id": object})
        intact_df = intact_df[~intact_df["pubmed_id"].isin(params.pid_to_drop)] # Remove merged bioplex

        ppi_cvcl_df = intact_df.merge(single_pid_df, on="pubmed_id", how="left")
        ppi_cvcl_df = ppi_cvcl_df.merge(manual_df, on="pubmed_id", how="left", suffixes=("_pubtator", "_manual"))
        get_cvcl = lambda x: x["CVCL_manual"] if x["CVCL_manual"] == x["CVCL_manual"] else x["CVCL_pubtator"] # remove nan and those single that i accidently checked
        ppi_cvcl_df["CVCL"] = ppi_cvcl_df.apply(get_cvcl, axis=1)
        ppi_cvcl_df = ppi_cvcl_df[~ppi_cvcl_df["CVCL"].isna()]

        # Get_bioplex
        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        bioplex_0291 = format_bioplex(input.cvcl_0291_bp, gene_name_df, "CVCL_0291")
        bioplex_0063 = format_bioplex(input.cvcl_0063_bp, gene_name_df, "CVCL_0063")

        all_cvcl_ppi = pd.concat([
            ppi_cvcl_df[column_order],
            bioplex_0291[column_order],
            bioplex_0063[column_order]
        ])
        all_cvcl_ppi.to_csv(output.cvcl_ppi, sep="\t", index=False)