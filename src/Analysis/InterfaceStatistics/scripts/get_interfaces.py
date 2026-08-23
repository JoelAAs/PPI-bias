import gzip
import os
from collections import defaultdict
import requests
from Bio.PDB import MMCIFParser, NeighborSearch
import pandas as pd


def build_entry_accessions(sifts_path, human_proteins, log):
    ent = defaultdict(set)
    ent_chains = defaultdict(lambda: defaultdict(set))
    opener = gzip.open if sifts_path.endswith(".gz") else open
    n = 0
    with opener(sifts_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3 or f[2] == "SP_PRIMARY":
                continue
            pdb, chain, acc = f[0].upper(), f[1], f[2]
            if acc not in human_proteins:
                continue

            ent[pdb].add(acc)
            ent_chains[pdb][chain].add(acc)
            n += 1

    print(f"Uniprot->PDB ids: {n} from {sifts_path}", file=log, flush=True)
    ent = {e: a for e, a in ent.items() if len(a) >= 2}
    ent_chains = {e: dict(ent_chains[e]) for e in ent}
    print(f"PDB entries with >=2 human chains: {len(ent)}", file=log, flush=True)

    return ent, ent_chains


def load_protein_pair_file(pod_file):
    df = pd.read_parquet(pod_file)
    proteins = set(df["uniprot_id_bait"]) | set(df["uniprot_id_prey"])
    return proteins, df


def match_protein_pairs_to_pdbs(pp_df, ent, session, log):
    prot2ent_dict = defaultdict(set)
    for ent_id, prot_ids in ent.items():
        for prot in prot_ids:
            prot2ent_dict[prot].add(ent_id)

    pp_df["matching_pdb"] = pp_df.apply(
        lambda row: prot2ent_dict[row["uniprot_id_bait"]]
        & prot2ent_dict[row["uniprot_id_prey"]],
        axis=1,
    )

    pp_df = pp_df[pp_df["matching_pdb"].map(len) > 0]
    print(f"Pairs with a matching PDB entry: {len(pp_df)}", file=log, flush=True)
    pp_df["matching_pdb"] = pp_df["matching_pdb"].map(
        lambda x: get_min_resolution(x, session)
    )
    print("Resolved best-resolution PDB entry per pair", file=log, flush=True)

    return pp_df


def get_min_resolution(pdb_ids, session):
    if len(pdb_ids) == 1:
        return list(pdb_ids)[0]

    min_res = 9000
    top_id = None
    for pdb_id in pdb_ids:
        resolution = get_resolution(pdb_id, session)
        if resolution < min_res:
            min_res = resolution
            top_id = pdb_id

    return top_id


def get_resolution(pdb_id, session):
    response = session.get(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=30
    )

    if response.status_code != 200:
        raise RuntimeError(
            f"Failed to fetch PDB metadata for {pdb_id}: status {response.status_code}"
        )

    pdb_content = response.json()
    resolution = pdb_content.get("rcsb_entry_info", {}).get("resolution_combined")
    return resolution


def download_pdb_structures(pp_df, session, pdb_folder, log):
    os.makedirs(pdb_folder, exist_ok=True)

    pdb_ids = pp_df["matching_pdb"].unique()
    print(f"Downloading {len(pdb_ids)} PDB structures", file=log, flush=True)
    for i, pdb_id in enumerate(pdb_ids, start=1):
        path = os.path.join(pdb_folder, f"{pdb_id}.cif.gz")
        response = session.get(
            f"https://files.rcsb.org/download/{pdb_id}.cif.gz", timeout=30
        )

        if response.status_code != 200:
            raise RuntimeError(
                f"Failed to download PDB structure for {pdb_id}: status {response.status_code}"
            )

        with open(path, "wb") as f:
            f.write(response.content)

        if i % 50 == 0 or i == len(pdb_ids):
            print(f"Downloaded {i}/{len(pdb_ids)} PDB structures", file=log, flush=True)


def get_interface_residues(structure, chains_a, chains_b, max_distance=6.0):
    model = next(structure.get_models())

    atoms_a = [
        atom for chain in model if chain.id in chains_a for atom in chain.get_atoms()
    ]
    atoms_b = [
        atom for chain in model if chain.id in chains_b for atom in chain.get_atoms()
    ]

    ns_b = NeighborSearch(atoms_b)
    ns_a = NeighborSearch(atoms_a)
    n_residues = 0

    res_b = set()
    for atom in atoms_a:
        for r in ns_b.search(atom.coord, max_distance, level="R"):
            res_b.add(r)

    res_a = set()
    for atom in atoms_b:
        for r in ns_a.search(atom.coord, max_distance, level="R"):
            res_a.add(r)

    interface_size = (len(res_a) + len(res_b)) / 2
    return interface_size


def get_interface_size(pdb_id, acc_a, acc_b, ent_chains, pdb_folder, max_distance=6.0):
    chains_a = [chain for chain, prots in ent_chains[pdb_id].items() if acc_a in prots]
    chains_b = [chain for chain, prots in ent_chains[pdb_id].items() if acc_b in prots]
    if not chains_a or not chains_b:
        raise ValueError(f"Chains for {acc_a} or {acc_b} not found in PDB {pdb_id}")

    with gzip.open(f"{pdb_folder}/{pdb_id}.cif.gz", "rt") as fh:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, fh)

    interface_size = get_interface_residues(structure, chains_a, chains_b, max_distance)
    return interface_size


def set_interface_sizes(pp_df, ent_chains, pdb_folder, log):
    pp_df["interface_residues"] = pp_df.apply(
        lambda row: get_interface_size(
            row["matching_pdb"],
            row["uniprot_id_bait"],
            row["uniprot_id_prey"],
            ent_chains,
            pdb_folder,
        ),
        axis=1,
    )
    print(f"Computed interface sizes for {len(pp_df)} pairs", file=log, flush=True)
    return pp_df


if __name__ == "__main__":
    pod_file = snakemake.input.pod_file
    sifts_path = snakemake.input.sifts_file
    pdb_folder = snakemake.output.pdb_folder
    output_file = snakemake.output.interface_annotated
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    human_proteins, pp_df = load_protein_pair_file(pod_file)
    print(f"Loaded {len(pp_df)} pairs over {len(human_proteins)} proteins", file=log, flush=True)
    ent, ent_chains = build_entry_accessions(sifts_path, human_proteins, log)

    with requests.Session() as session:
        pp_df = match_protein_pairs_to_pdbs(pp_df, ent, session, log)
        download_pdb_structures(pp_df, session, pdb_folder, log)

    pp_df = set_interface_sizes(pp_df, ent_chains, pdb_folder, log)
    pp_df.to_parquet(output_file, index=False)
    print("Done.", file=log, flush=True)
    log.close()
