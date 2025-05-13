import pandas as pd
import re
import numpy as np

## TODO: separate and rename files, this isn't all formating


def _find_pattern(cell, pattern, single=True):
    match = re.search(pattern, cell)
    if match:
        matches = match.groups()
        if single:
            return matches[0]
        return match.groups()
    return None

def get_gene_names(filename, output_file):
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

def read_intact(filename):
    intact_df = pd.read_csv(filename, sep="\t")
    intact_df = intact_df[
        (
                intact_df['Taxid interactor A'].str.contains("9606") &
                intact_df['Taxid interactor B'].str.contains("9606")
        )
    ]

    intact_df['RoleA'] = intact_df['Experimental role(s) interactor A'].str.contains("bait")
    intact_df['RoleB'] = intact_df['Experimental role(s) interactor B'].str.contains("bait")
    # No bait-bait prey-prey pairs
    intact_df = intact_df[intact_df['RoleA'] ^ intact_df['RoleB']]


    intact_df["pubmedID"] = intact_df["Publication Identifier(s)"].apply(
        _find_pattern, args=(r"pubmed:(\d+)",))
    intact_df["IDA"] = intact_df["#ID(s) interactor A"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    intact_df["IDB"] = intact_df["ID(s) interactor B"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    intact_df["detectionMethod"] = intact_df["Interaction detection method(s)"].apply(
        _find_pattern, args=(r"psi-mi:\"(MI:\d+)\"",))


    return intact_df

def filter_and_index(intact_df):
    intact_df = intact_df[["IDA", "IDB", "RoleA", "RoleB", "detectionMethod", "pubmedID"]].drop_duplicates(keep="first")
    intact_df = intact_df.dropna()
    intact_df = intact_df[intact_df["RoleA"] ^ intact_df["RoleB"]]
    reform_list = []
    for _, row in intact_df.iterrows():
        if row["RoleA"]:
            bait = row["IDA"]
            prey = row["IDB"]
        else:
            bait = row["IDB"]
            prey = row["IDA"]

        reform_list.append({
            "bait": bait,
            "prey": prey,
            "pubmed_id": row["pubmedID"],
            "detection_method": row["detectionMethod"].replace(":", "-")
        })

    bait_prey_df = pd.DataFrame(reform_list)
    bait_prey_df = bait_prey_df.drop_duplicates()

    return bait_prey_df

def create_or_update(current_dict, bait, prey):
    if bait in current_dict:
        if prey in current_dict[bait]:
            current_dict[bait][prey] += 1
        else:
            current_dict[bait][prey] = 1
    else:
        current_dict[bait] = {prey: 1}

    return current_dict


def get_interaction_dict(bait_prey_df, method="Y2H-pooling", prey_file=None):
    cell_line_preys = []
    if method == "MS":
        prey_df = pd.read_csv(prey_file, sep="\t")
        cell_line_preys = prey_df["uniprot_id"].tolist()
    elif method != "Y2H-pooling":
        raise ValueError(f"The method type: {method} is has no prey-bait matching strategy")


    max_interaction_dict = dict()
    observed_interaction_dict = dict()

    studies = bait_prey_df["pubmed_id"].unique()
    for study in studies:
        study_df = bait_prey_df[bait_prey_df["pubmed_id"] == study]
        baits = set()
        preys = set()
        for _, row in study_df.iterrows():
            bait = row["bait"]
            prey = row["prey"]

            baits.update({bait})
            preys.update({prey})

            observed_interaction_dict = create_or_update(
                observed_interaction_dict,
                bait=bait,
                prey=prey
            )

            if method == "MS" and prey not in cell_line_preys:
                # Identified prey without POD estimate
                max_interaction_dict = create_or_update(
                    max_interaction_dict,
                    bait=bait,
                    prey=prey
                )
        if method == "MS":
            # Observe we assume that test all baits as separate MS runs.
            # This should be the standard way of doing things.
            preys = cell_line_preys

        for bait in baits:
            for prey in preys:
                max_interaction_dict = create_or_update(
                    max_interaction_dict,
                    bait=bait,
                    prey=prey
                )

    return max_interaction_dict, observed_interaction_dict

def dict_to_pairs_file(max_interaction_dict, observed_interaction_dict, output_filename, method = "Y2H-pooling", prey_pod_file =None):
    pod_dict = dict()
    if method == "MS":
        with open(prey_pod_file, "r") as f:
            for l in f:
                gene_name, protein_id, prey_pod = l.strip().split("\t")
                pod_dict[protein_id] = prey_pod


    with open(output_filename, "w") as w:
        w.write("prey\tbait\tmax_interactions\tobserved_interactions\tprey_pod\n")
        for bait in max_interaction_dict:
            for prey in max_interaction_dict[bait]:
                if method == "MS":
                    if prey in pod_dict:
                        prey_pod = pod_dict[prey]
                    else:
                        prey_pod = None
                else:
                    prey_pod = 1

                obs_count = 0
                if prey in observed_interaction_dict[bait]:
                    obs_count = observed_interaction_dict[bait][prey]

                w.write(f"{bait}\t"
                        f"{prey}\t"
                        f"{max_interaction_dict[bait][prey]}\t"
                        f"{obs_count}\t"
                        f"{prey_pod}\n")


