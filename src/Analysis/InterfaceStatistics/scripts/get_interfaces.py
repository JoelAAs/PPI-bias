import gzip
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, NeighborSearch
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from datetime import datetime


def make_session(total_retries=5, backoff_factor=1.0):
    retry = Retry(
        total=total_retries,
        connect=total_retries,
        read=total_retries,
        status=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


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
    pp_df = pp_df[pp_df["matching_pdb"].map(len) > 0].copy()
    print(f"Pairs with a matching PDB entry: {len(pp_df)}", file=log, flush=True)
    unique_multi_pdbs = set().union(*[pdbs for pdbs in pp_df["matching_pdb"] if len(pdbs) > 1])

    resolution_dict = {pdb: get_resolution(pdb, session) for pdb in unique_multi_pdbs}

    print("Resolved best-resolution PDB entry per pair", file=log, flush=True)

    def _get_min_resolution(pdb_set, resolution_dict):
        if len(pdb_set) == 1:
            return next(iter(pdb_set))

        min_res = float("inf")
        best_pdb = None
        for pdb in pdb_set:
            res = resolution_dict.get(pdb)
            if res is None: #nmr
                res = 50 # set large
            if res < min_res:
                min_res = res
                best_pdb = pdb
        return best_pdb

    pp_df["matching_pdb"] = pp_df["matching_pdb"].apply(lambda x: _get_min_resolution(x, resolution_dict))

    return pp_df


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
    if not resolution:
        return None
    return min(resolution)


def download_pdb_structures(pp_df, session, pdb_folder, log):
    os.makedirs(pdb_folder, exist_ok=True)

    pdb_ids = pp_df["matching_pdb"].unique()
    print(f"Downloading {len(pdb_ids)} PDB structures", file=log, flush=True)
    n_cached = 0
    for i, pdb_id in enumerate(pdb_ids, start=1):
        path = os.path.join(pdb_folder, f"{pdb_id}.cif.gz")
        if os.path.exists(path):
            n_cached += 1
            continue

        s = datetime.now()
        response = session.get(
            f"https://files.rcsb.org/download/{pdb_id}.cif.gz", timeout=30
        )

        if response.status_code != 200:
            print(
                f"Failed to download PDB structure for {pdb_id}: status {response.status_code}",
                file=log, flush=True,
            )
            continue

        with open(path, "wb") as f:
            f.write(response.content)
        e = datetime.now()
        print(f"Downloaded {pdb_id} in {(e - s).total_seconds():.2f} seconds", file=log, flush=True)

        if i % 50 == 0 or i == len(pdb_ids):
            print(f"Progress: {i}/{len(pdb_ids)} PDB structures ({n_cached} already cached)", file=log, flush=True)


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
            
    n_res_a, n_res_b = len(res_a), len(res_b)
    avg = (n_res_a + n_res_b) / 2
    return n_res_a, n_res_b, avg

def get_interface_sizes_for_structure(pdb_id, pairs, ent_chains, pdb_folder, min_residues_per_side, max_distance=6.0):
    with gzip.open(f"{pdb_folder}/{pdb_id}.cif.gz", "rt") as fh:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, fh)

    sizes = {}
    for idx, acc_a, acc_b in pairs:
        chains_a = [chain for chain, prots in ent_chains[pdb_id].items() if acc_a in prots]
        chains_b = [chain for chain, prots in ent_chains[pdb_id].items() if acc_b in prots]
        if not chains_a or not chains_b:
            raise ValueError(f"Chains for {acc_a} or {acc_b} not found in PDB {pdb_id}")
        n_res_a, n_res_b, avg = get_interface_residues(structure, chains_a, chains_b, max_distance)
        if min(n_res_a, n_res_b) < min_residues_per_side:
            sizes[idx] = 0
        else:
            sizes[idx] = avg
    return sizes


def set_interface_sizes(pp_df, ent_chains, pdb_folder, log, min_residues_per_side=3, max_workers=None):
    # check residue contanct on both sides of the interface
    groups = [
        (pdb_id, list(sub_df[["uniprot_id_bait", "uniprot_id_prey"]].itertuples(name=None)))
        for pdb_id, sub_df in pp_df.groupby("matching_pdb")
    ]

    sizes = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(get_interface_sizes_for_structure, pdb_id, pairs, ent_chains, pdb_folder, min_residues_per_side)
            for pdb_id, pairs in groups
        ]
        for future in futures:
            sizes.update(future.result())

    pp_df["interface_residues"] = pp_df.index.map(sizes)
    pp_df = pp_df[pp_df["interface_residues"] > 0] # remoce pairs with no interface residues on either side
    print(f"Computed interface sizes for {len(pp_df)} pairs across {len(groups)} structures", file=log, flush=True)
    return pp_df


if __name__ == "__main__":
    matched_pairs_file = snakemake.input.matched_pairs
    sifts_path = snakemake.input.sifts_file
    pdb_folder = snakemake.input.pdb_folder
    output_file = snakemake.output.interface_annotated
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    pp_df = pd.read_parquet(matched_pairs_file)
    human_proteins = set(pp_df["uniprot_id_bait"]) | set(pp_df["uniprot_id_prey"])
    print(f"Loaded {len(pp_df)} matched pairs over {len(human_proteins)} proteins", file=log, flush=True)
    _, ent_chains = build_entry_accessions(sifts_path, human_proteins, log)

    pp_df = set_interface_sizes(pp_df, ent_chains, pdb_folder, log, max_workers=snakemake.threads)
    pp_df.to_parquet(output_file, index=False)
    print("Done.", file=log, flush=True)
    log.close()
