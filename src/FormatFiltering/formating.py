import gzip
import re
import time
import warnings

import pandas as pd
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


def parse_sec_ac_mapping(sec_ac_file):
    """
    Parse UniProt sec_ac.txt into {secondary_ac: primary_ac}.
    Identifies data lines by requiring both tokens to match the UniProt accession pattern.
    """
    ac_pattern = re.compile(r'^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$')
    mapping = {}
    with open(sec_ac_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2 and ac_pattern.match(parts[0]) and ac_pattern.match(parts[1]):
                mapping[parts[0]] = parts[1]
    return mapping


def _get_displayed_isoforms(base_acs):
    """
    Query UniProt for the displayed (canonical) isoform accession of each base AC.
    Returns {base_ac: displayed_isoform_accession} for entries that have explicit
    alternative products. Base ACs absent from the result have no isoform annotation
    and are canonical as bare accessions.
    """
    displayed = {}
    base_acs = list(base_acs)
    for i in range(0, len(base_acs), 100):
        batch = base_acs[i:i + 100]
        query = " OR ".join(f"accession:{ac}" for ac in batch)
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": query,
                "format": "json",
                "fields": "accession,cc_alternative_products",
                "size": len(batch) + 10,
            },
        )
        if not resp.ok:
            continue
        for entry in resp.json().get("results", []):
            ac = entry["primaryAccession"]
            for comment in entry.get("comments", []):
                if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
                    for isoform in comment.get("isoforms", []):
                        if isoform.get("isoformSequenceStatus") == "Displayed":
                            ids = isoform.get("isoformIds", [])
                            if ids:
                                displayed[ac] = ids[0]
                            break
                    break
        time.sleep(0.05)
    return displayed


def build_entrez_mapping(mitab_file, sec_ac_file, idmapping_gz, gene_info_gz, output_file, unmapped_file=None, keep_non_canonical=False):
    """
    Map canonical UniProt IDs in miTab to numeric Entrez Gene IDs.
    Isoform-suffixed IDs are kept only when UniProt marks them as the displayed
    (canonical) isoform; others are excluded so those interactions are dropped
    by the inner join in format_miTab.

    Steps:
    1. Collect all UniProt IDs from the human-only rows of miTab.
    2. For isoform-suffixed IDs, query UniProt to find the displayed isoform;
       keep only those that are canonical, drop the rest.
    3. Remap secondary/retired ACs to current primary ACs via sec_ac.txt.
    4. Look up Entrez GeneID in HUMAN_9606_idmapping.dat.gz.
    5. For entries missing a direct GeneID: resolve via gene symbol in
       Homo_sapiens.gene_info.gz (Symbol and Synonyms columns).
    """
    intact_df = pd.read_csv(mitab_file, sep="\t")
    intact_df = intact_df[
        intact_df['Taxid interactor A'].str.contains("9606") &
        intact_df['Taxid interactor B'].str.contains("9606")
    ]
    intact_df["IDA"] = intact_df["#ID(s) interactor A"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))
    intact_df["IDB"] = intact_df["ID(s) interactor B"].apply(
        _find_pattern, args=(r"uniprotkb:(.+)",))

    raw_ids = {i for i in (list(intact_df["IDA"]) + list(intact_df["IDB"])) if i}

    bare_ids = {uid for uid in raw_ids if "-" not in uid}
    isoform_ids = {uid for uid in raw_ids if "-" in uid}

    canonical_isoforms = set()
    if isoform_ids:
        base_acs_with_isoforms = {uid.split("-")[0] for uid in isoform_ids}
        displayed = _get_displayed_isoforms(base_acs_with_isoforms)
        # Keep an isoform ID only when it is the displayed isoform for its base AC
        for uid in isoform_ids:
            base = uid.split("-")[0]
            if displayed.get(base) == uid or keep_non_canonical:
                canonical_isoforms.add(uid)

    dropped = len(isoform_ids) - len(canonical_isoforms)
    if dropped and not keep_non_canonical:
        warnings.warn(f"{dropped} non-canonical isoforms excluded (interactions dropped).")

    all_ids = list(bare_ids | canonical_isoforms)

    sec_ac_map = parse_sec_ac_mapping(sec_ac_file)

    id_to_primary = {}
    for uid in all_ids:
        base = uid.split("-")[0]
        id_to_primary[uid] = sec_ac_map.get(base, base)

    lookup_ids = set(id_to_primary.values())

    uniprot_to_entrez = {}
    uniprot_to_symbol = {}
    with gzip.open(idmapping_gz, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            ac, id_type, value = parts
            if ac not in lookup_ids:
                continue
            if id_type == "GeneID":
                uniprot_to_entrez[ac] = value
            elif id_type == "Gene_Name" and ac not in uniprot_to_symbol:
                uniprot_to_symbol[ac] = value

    missing_acs = [ac for ac in lookup_ids if ac not in uniprot_to_entrez]
    if missing_acs:
        missing_symbols = {uniprot_to_symbol[ac] for ac in missing_acs if ac in uniprot_to_symbol}
        symbol_to_entrez = {}
        with gzip.open(gene_info_gz, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                gene_id, symbol, synonyms_field = parts[1], parts[2], parts[4]
                synonyms = synonyms_field.split("|") if synonyms_field != "-" else []
                if symbol in missing_symbols:
                    symbol_to_entrez[symbol] = gene_id
                for syn in synonyms:
                    if syn in missing_symbols:
                        symbol_to_entrez.setdefault(syn, gene_id)

        for ac in missing_acs:
            sym = uniprot_to_symbol.get(ac)
            if sym and sym in symbol_to_entrez:
                uniprot_to_entrez[ac] = symbol_to_entrez[sym]

    rows = []
    dropped = []
    for uid in all_ids:
        primary = id_to_primary[uid]
        entrez = uniprot_to_entrez.get(primary)
        if entrez:
            rows.append({"uniprot_id": uid, "gene_name": entrez})
        else:
            dropped.append({"uniprot_id": uid, "primary_ac": primary})

    if dropped:
        warnings.warn(f"{len(dropped)} UniProt IDs had no Entrez GeneID mapping and were dropped.")

    pd.DataFrame(rows).drop_duplicates().to_csv(output_file, sep="\t", index=False)

    if unmapped_file is not None:
        pd.DataFrame(dropped).to_csv(unmapped_file, sep="\t", index=False)
