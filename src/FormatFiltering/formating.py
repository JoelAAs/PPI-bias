import time
import pandas as pd
import re
import requests


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
            "uniprot_id_bait": bait,
            "uniprot_id_prey": prey,
            "pubmed_id": row["pubmed_id"],
            "detection_method": row["detection_method"].replace(":", "-")
        })

    bait_prey_df = pd.DataFrame(reform_list)
    bait_prey_df = bait_prey_df.drop_duplicates()

    return bait_prey_df


def get_gene_names(filename, output_file):
    """
    Primary gene name registered for accession id in uniprot, those without are discarded
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


    all_ids = list(set(intact_df["IDB"]) | set(intact_df["IDA"]))
    all_ids = [id for id in all_ids if id]

    def binit(x, n):
        for i in range(0,len(x), n):
            yield x[i:i+n]

    id_bins = binit(all_ids, 100)

    with open(output_file, "w") as w:
        url = "https://rest.uniprot.org/uniprotkb/search"

        w.write("uniprot_id\tgene_name\n")
        for bin in id_bins:
            params = {
                "query": f"accession:{' OR '.join(bin)}",
                "format": "tsv",
                "fields": "accession,gene_primary",
                "size":200
            }
            response = requests.get(url, params=params)
            if response.ok:
                lines = response.text.split("\n")

                for line in lines[1:]:
                    if line:
                        w.write(line + "\n")

            time.sleep(0.01)


    gene_name_df = pd.read_csv(output_file, sep="\t")
    gene_name_df = gene_name_df.drop_duplicates()
    missing_queries = [a for a in all_ids if a not in gene_name_df.uniprot_id.values]
    isoform_rows = []
    for missing in missing_queries:
        if "-" in missing:
            p_id, isoform_n = missing.split("-")
            gene = gene_name_df[gene_name_df["uniprot_id"] == p_id]["gene_name"]
            if not gene.empty:
                gene = gene.values[0]
                isoform_rows.append({
                    "uniprot_id": missing,
                    "gene_name": f"{gene}-{isoform_n}"
                })

    combined_df = pd.concat([
        gene_name_df,
        pd.DataFrame(isoform_rows)
    ], ignore_index=True)
    combined_df.to_csv(output_file, sep="\t", index=False)

