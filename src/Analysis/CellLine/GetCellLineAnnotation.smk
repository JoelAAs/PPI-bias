import pandas as pd
from src.FormatFiltering.formating import parse_sec_ac_mapping

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


rule get_pubtator:
    output:
        "work_folder/cell_line_annotation/cellline2pubtator3.txt"
    log:
        "logs/cell_line_annotation/cellline2pubtator3.log"
    shell:
        """
        exec > {log} 2>&1
        wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/cellline2pubtator3.gz
        gunzip cellline2pubtator3.gz
        mv cellline2pubtator3 {output}
        """


rule get_bioplex:
    output:
        cvcl_0291_bp="work_folder/data/bioplex/CVCL_0291_bp.csv",
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063_bp.csv",
        cvcl_0291_net="work_folder/data/bioplex/CVCL_0291_net.csv",
        cvcl_0063_net="work_folder/data/bioplex/CVCL_0063_net.csv",
    log:
        "logs/data/bioplex/bioplex.log"
    shell:
        """
        exec > {log} 2>&1
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_293T_DirectedEdges.tsv -O {output.cvcl_0063_bp}
        wget https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv -O {output.cvcl_0063_net}
        wget https://bioplex.hms.harvard.edu/data/BioPlex_3.0_HCT116_DirectedEdges.tsv -O {output.cvcl_0291_bp}
        wget https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv -O {output.cvcl_0291_net}
        """



rule format_bioplex:
    input:
        bp="work_folder/data/bioplex/{bp_cl}_bp.csv",
        network="work_folder/data/bioplex/{bp_cl}_net.csv",
        sec_ac="work_folder/data/uniprot/sec_ac.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        formated="work_folder/data/bioplex/{bp_cl}.csv",
    log:
        "logs/data/bioplex/{bp_cl}_format.log"
    run:
        df_bp = pd.read_csv(input.bp, sep="\t")[["Bait GeneID", "Prey GeneID"]]
        df_bp.columns = ["gene_id_bait", "gene_id_prey"]
        df_network = pd.read_csv(input.network, sep="\t")[["GeneA", "GeneB", "UniprotA", "UniprotB"]]

        net_fwd = df_network.rename(columns={
            "GeneA": "gene_id_bait", "GeneB": "gene_id_prey",
            "UniprotA": "uniprot_id_bait", "UniprotB": "uniprot_id_prey"})
        net_rev = df_network.rename(columns={
            "GeneB": "gene_id_bait", "GeneA": "gene_id_prey",
            "UniprotB": "uniprot_id_bait", "UniprotA": "uniprot_id_prey"})
        net_both = pd.concat([net_fwd, net_rev]).drop_duplicates(["gene_id_bait", "gene_id_prey"])

        df_joined = df_bp.merge(net_both, on=["gene_id_bait", "gene_id_prey"], how="left", validate="many_to_one")

        n_missing = df_joined["uniprot_id_bait"].isna().sum()
        assert n_missing == 0, f"{n_missing} bp edges not found in network"
        assert len(df_joined) == len(df_bp), "Join changed number of bp edges"

        net_keys = set(tuple(sorted(s)) for s in zip(df_network["GeneA"], df_network["GeneB"]))
        bp_keys = set(tuple(sorted(s)) for s in zip(df_bp["gene_id_bait"], df_bp["gene_id_prey"]))
        
        assert net_keys == bp_keys, f"{len(net_keys - bp_keys)} network edges not in bp, {len(bp_keys - net_keys)} bp edges not in network"

        ac_map = parse_sec_ac_mapping(input.sec_ac)
        
        def _map_ac(ac_id):
            iso_n = None
            if "-" in ac_id:
                ac_id_s, iso_n = ac_id.split("-")
            else:
                ac_id_s = ac_id

            if ac_id_s in ac_map:
                if iso_n:
                    return f"{ac_map[ac_id_s]}-{iso_n}"
                return ac_map[ac_id_s]
            return ac_id

        df_joined["uniprot_id_bait"] = df_joined["uniprot_id_bait"].apply(_map_ac)
        df_joined["uniprot_id_prey"] = df_joined["uniprot_id_prey"].apply(_map_ac)


        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        df_joined = df_joined.merge(gene_name_df, left_on="uniprot_id_bait", right_on="uniprot_id")
        del df_joined["uniprot_id"]
        df_joined = df_joined.merge(gene_name_df, left_on="uniprot_id_prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
        del df_joined["uniprot_id"]

        df_joined["pubmed_id"] = "33961781"
        df_joined["detection_method"] = "MI-0007"
        df_joined["CVCL"] = wildcards.bp_cl
        column_order =[
            'uniprot_id_bait', 'uniprot_id_prey', 'pubmed_id', 'detection_method',
            'gene_name_bait', 'gene_name_prey', 'CVCL'
        ]
        df_joined[column_order].to_csv(output.formated, sep="\t", index=False)
        


rule get_cellosaurus:
    output:
        cellosaurus = "work_folder/data/cellosaurus/cellosaurus.txt"
    log:
        "logs/data/cellosaurus/cellosaurus.log"
    shell:
        """
        wget https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt -O {output.cellosaurus} > {log} 2>&1
        """

rule get_cellosaurus_human_cl:
    input:
        cellosaurus = "work_folder/data/cellosaurus/cellosaurus.txt"
    output:
        cellosaurus_csv = "work_folder/data/cellosaurus/taxon_cellosaurus.csv",
        cellosaurus_human = "work_folder/data/cellosaurus/cellosaurus_human.txt"
    log:
        "logs/data/cellosaurus/cellosaurus_human.log"
    run:
        parse_cellosaurus(input.cellosaurus, output.cellosaurus_csv)
        df_cellosaurus = pd.read_csv(output.cellosaurus_csv, sep="\t")
        df_cellosaurus = df_cellosaurus[df_cellosaurus["OX"].str.contains("9606")]
        df_cellosaurus["AC"].to_csv(output.cellosaurus_human, index=None)




rule join_to_pid:
    params:
        pid_to_drop = config["pid_to_remove_cl"]
    input:
        pubator = "work_folder/cell_line_annotation/cellline2pubtator3.txt",
        human_cell_lines = "work_folder/data/cellosaurus/cellosaurus_human.txt",
        formatted_intact = "work_folder/formated/bait_prey_publications.csv",
        cvcl_0291_bp="work_folder/data/bioplex/CVCL_0291.csv",
        cvcl_0063_bp="work_folder/data/bioplex/CVCL_0063.csv",
        manual_curated = "data/checked_studies.csv"
    output:
        cvcl_ppi = "work_folder/formated/bait_prey_CVCL.csv"
    log:
        "logs/formated/bait_prey_CVCL.log"
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
        bioplex_0291 = pd.read_csv(input.cvcl_0291_bp, sep="\t")
        bioplex_0063 = pd.read_csv(input.cvcl_0063_bp, sep="\t")

        all_cvcl_ppi = pd.concat([
            ppi_cvcl_df[column_order],
            bioplex_0291[column_order],
            bioplex_0063[column_order]
        ])
        all_cvcl_ppi.to_csv(output.cvcl_ppi, sep="\t", index=False)