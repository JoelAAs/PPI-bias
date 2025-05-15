import pandas as pd
import re
import numpy as np

def _find_pattern(cell, pattern, single=True):
    match = re.search(pattern, cell)
    if match:
        matches = match.groups()
        if single:
            return matches[0]
        return match.groups()
    return None

def filter_mitab(filename):
    """
    Remove non-human proteins.
    Remove any non bait-prey interaction.
    Format publications ids and uniprot identifiers
    :param filename: miTab interaction file
    :return: filtered and formated pandas dataframe
    """
    mitab_df = pd.read_csv(filename, sep="\t")
    mitab_df = mitab_df[
        (
                mitab_df['Taxid interactor A'].str.contains("9606") &
                mitab_df['Taxid interactor B'].str.contains("9606")
        )
    ]
    mitab_df['RoleA'] = mitab_df['Experimental role(s) interactor A'].str.contains("bait")
    mitab_df['RoleB'] = mitab_df['Experimental role(s) interactor B'].str.contains("bait")

    # No bait-bait prey-prey pairs
    mitab_df = mitab_df[mitab_df['RoleA'] ^ mitab_df['RoleB']]

    mitab_df["pubmed_id"] = mitab_df["Publication Identifier(s)"].apply(
        _find_pattern, args=(r"pubmed:(\d+)",))
    mitab_df["IDA"] = mitab_df["#ID(s) interactor A"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    mitab_df["IDB"] = mitab_df["ID(s) interactor B"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    mitab_df["detection_method"] = mitab_df["Interaction detection method(s)"].apply(
        _find_pattern, args=(r"psi-mi:\"(MI:\d+)\"",))

    return mitab_df


def reform_to_bait_prey(mitab_df):
    """
    Reform PPI pandas dataframe into bait-prey format
    :param mitab_df: input filtered
    :return: data frame with columns:
        uniprot_bait
        uniprot_prey
        pubmed_id
        detection_method
    """
    mitab_df = mitab_df[["IDA", "IDB", "RoleA", "RoleB", "detection_method", "pubmed_id"]].drop_duplicates(keep="first")
    mitab_df = mitab_df.dropna()
    mitab_df = mitab_df[mitab_df["RoleA"] ^ mitab_df["RoleB"]]
    reform_list = []
    for _, row in mitab_df.iterrows():
        if row["RoleA"]:
            bait = row["IDA"]
            prey = row["IDB"]
        else:
            bait = row["IDB"]
            prey = row["IDA"]

        reform_list.append({
            "uniprot_bait": bait,
            "uniprot_prey": prey,
            "pubmed_id": row["pubmed_id"],
            "detection_method": row["detection_method"].replace(":", "-")
        })

    bait_prey_df = pd.DataFrame(reform_list)
    bait_prey_df = bait_prey_df.drop_duplicates()

    return bait_prey_df


def get_gene_names(filename, output_file):
    """
    Extracting uniprot-gene name from miTab interaction file
    :param filename: Mitab file location
    :param output_file: uniprot - gene_name file location
    :return:
    """
    intact_df = pd.read_csv(filename, sep="\t")
    intact_df = intact_df[
        (
                intact_df['Taxid interactor A'].str.contains("9606") &
                intact_df['Taxid interactor B'].str.contains("9606")
        )
    ]
    intact_df["IDA"] = intact_df["#ID(s) interactor A"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    intact_df["IDB"] = intact_df["ID(s) interactor B"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))

    gene_name_pattern = r"uniprotkb:([a-zA-Z0-9-]+)\(gene name\)"
    intact_df["gene_name_a"] = intact_df["Alias(es) interactor A"].apply(
        _find_pattern, args=(gene_name_pattern,))
    intact_df["gene_name_b"] = intact_df["Alias(es) interactor B"].apply(
        _find_pattern, args=(gene_name_pattern,))

    name_dict = dict()
    for i, row in intact_df.iterrows():
        for uniprot_id, gene_name in zip(
                [row["IDA"], row["IDB"]],
                [row["gene_name_a"], row["gene_name_b"]]
        ):
            if uniprot_id and gene_name and uniprot_id not in name_dict:
                name_dict[uniprot_id] = gene_name

    with open(output_file, "w") as w:
        w.write("uniprot_id\tgene_name\n")
        for uniprot_id, gene_name in name_dict.items():
            w.write(f"{uniprot_id}\t{gene_name}\n")

