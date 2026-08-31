import gzip
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa
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


def build_chain_map(sifts_path, human_proteins, log):
    chain_accs = defaultdict(set)                            # (pdb,chain) -> {acc, ...}
    ranges = defaultdict(lambda: defaultdict(list))          # (pdb,chain) -> acc -> [(pdb_beg, pdb_end), ...]
    opener = gzip.open if sifts_path.endswith(".gz") else open
    n = 0
    with opener(sifts_path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\r\n").split("\t")
            if len(f) < 7 or f[2] == "SP_PRIMARY":
                continue
            acc = f[2]
            if acc not in human_proteins:
                continue
            pdb, chain = f[0].upper(), f[1]
            chain_accs[(pdb, chain)].add(acc)
            ranges[(pdb, chain)][acc].append((f[5], f[6]))
            n += 1

    def _residue_window(raw_ranges):
        allowed = set()
        for beg, end in raw_ranges:
            if not beg or not end:                          # blank author numbering, can't map
                continue
            try:
                b, e = int(beg), int(end)
            except ValueError:
                continue
            allowed.update(range(b, e + 1))
        return allowed

    chain_map = defaultdict(lambda: defaultdict(dict))
    fusions = []
    dropped = 0
    for (pdb, chain), accs in chain_accs.items():
        if len(accs) > 1:                                   # fusion: multiple accessions on one chain only use the largest segment
            windows = {a: _residue_window(ranges[(pdb, chain)][a]) for a in accs}
            dominant = max(windows, key=lambda a: len(windows[a]))
            fusions.append((pdb, chain, {a: len(w) for a, w in windows.items()}, dominant))
            if windows[dominant]:
                chain_map[pdb][dominant][chain] = windows[dominant]
            else:
                dropped += 1                                # no candidate has a usable residue range
        else:
            (acc,) = accs
            chain_map[pdb][acc][chain] = None

    chain_map = {p: {a: dict(c) for a, c in accs.items()} for p, accs in chain_map.items()}
    print(f"Uniprot->PDB ids: {n} from {sifts_path}", file=log, flush=True)
    print(f"Fusion chains detected: {len(fusions)} ({dropped} dropped, no usable residue window)", file=log, flush=True)
    print(f"PDB entries retained: {len(chain_map)}", file=log, flush=True)
    return chain_map, fusions


def read_fasta_lengths(fasta_file):
    lengths = {}
    with open(fasta_file) as f:
        acc = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    lengths[acc] = length
                acc = line.split("|")[1]
                length = 0
            else:
                length += len(line.strip())
        if acc is not None:
            lengths[acc] = length
    return lengths


def load_protein_pair_file(pod_file):
    df = pd.read_parquet(pod_file)
    proteins = set(df["uniprot_id_bait"]) | set(df["uniprot_id_prey"])
    return proteins, df


def match_protein_pairs_to_pdbs(pp_df, chain_map, session, log):
    """Pick one PDB entry per pair from chain_map (see build_chain_map).

    A candidate is only considered if both bait and prey are the kept
    (dominant) protein on at least one chain in that entry - a candidate
    where one side is only ever a fusion tag can never yield an interface,
    so it's excluded here rather than being selected and failing later.

    Among remaining candidates, entries are ranked first by how fusion-free
    they are for this specific pair (both sides on a clean, non-fusion
    chain beats one side fusion-restricted beats both sides restricted),
    and only resolution breaks ties within a tier - a worse-resolution
    clean structure is preferred over a better-resolution fusion structure.
    """
    prot2pdbs = defaultdict(set)
    for pdb, accs in chain_map.items():
        for acc in accs:
            prot2pdbs[acc].add(pdb)

    pp_df["matching_pdb"] = pp_df.apply(
        lambda row: prot2pdbs[row["uniprot_id_bait"]]
        & prot2pdbs[row["uniprot_id_prey"]],
        axis=1,
    )
    pp_df = pp_df[pp_df["matching_pdb"].map(len) > 0].copy()
    print(f"Pairs with a matching PDB entry: {len(pp_df)}", file=log, flush=True)
    unique_multi_pdbs = set().union(*[pdbs for pdbs in pp_df["matching_pdb"] if len(pdbs) > 1])

    resolution_dict = {pdb: get_resolution(pdb, session) for pdb in unique_multi_pdbs}

    print("Resolved best structure per pair (non-fusion first, then resolution)", file=log, flush=True)

    def _fusion_tier(pdb, acc):
        # 0 = has at least one clean (non-fusion) chain, 1 = fusion-restricted only
        windows = chain_map.get(pdb, {}).get(acc, {})
        return 0 if any(w is None for w in windows.values()) else 1

    def _select_pdb(row):
        pdb_set = row["matching_pdb"]
        if len(pdb_set) == 1:
            return next(iter(pdb_set))

        def _key(pdb):
            fusion_score = _fusion_tier(pdb, row["uniprot_id_bait"]) + _fusion_tier(pdb, row["uniprot_id_prey"])
            res = resolution_dict.get(pdb)
            if res is None: #nmr
                res = 50 # set large
            return (fusion_score, res)

        return min(pdb_set, key=_key)

    pp_df["matching_pdb"] = pp_df.apply(_select_pdb, axis=1)

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



def get_interface_residues(structure, smaller_windows, other_windows, max_distance=6.0):
    """smaller_windows/other_windows: {chain_id: allowed_resseq_set_or_None}.

    None means no residue-level restriction (use the whole chain); a set
    restricts atoms to that chain's real-protein residue range (fusion
    chains only, see build_chain_map).
    """
    model = next(structure.get_models())

    def _chain_atoms(windows):
        atoms = []
        for c in model:
            if c.id not in windows:
                continue
            allowed = windows[c.id]
            for res in c:
                if not is_aa(res, standard=False):
                    continue
                if allowed is not None and res.id[1] not in allowed:
                    continue
                atoms.extend(res.get_atoms())
        return atoms

    atoms_smaller = _chain_atoms(smaller_windows)
    atoms_other = _chain_atoms(other_windows)
    ns = NeighborSearch(atoms_smaller)
    positions = set()
    for atom in atoms_other:
        for r in ns.search(atom.coord, max_distance, level="R"):
            positions.add(r.id[1:])
    return len(positions)

def get_interface_sizes_for_structure(pdb_id, pairs, chain_map, pdb_folder, lengths, min_residues, max_distance=6.0):
    with gzip.open(f"{pdb_folder}/{pdb_id}.cif.gz", "rt") as fh:
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, fh)

    sizes = {}
    for idx, acc_a, acc_b in pairs:
        windows_a = chain_map.get(pdb_id, {}).get(acc_a, {})
        windows_b = chain_map.get(pdb_id, {}).get(acc_b, {})
        if not windows_a or not windows_b:
            sizes[idx] = 0
            continue

        if lengths[acc_a] <= lengths[acc_b]:
            smaller_windows_all, other_windows = windows_a, windows_b
        else:
            smaller_windows_all, other_windows = windows_b, windows_a

        n_residues = max(
            get_interface_residues(structure, {chain: allowed}, other_windows, max_distance)
            for chain, allowed in smaller_windows_all.items()
        )
        sizes[idx] = n_residues if n_residues >= min_residues else 0
    return sizes


def set_interface_sizes(pp_df, chain_map, pdb_folder, lengths, log, min_residues=3, max_workers=None):
    groups = [
        (pdb_id, list(sub_df[["uniprot_id_bait", "uniprot_id_prey"]].itertuples(name=None)))
        for pdb_id, sub_df in pp_df.groupby("matching_pdb")
    ]

    sizes = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(get_interface_sizes_for_structure, pdb_id, pairs, chain_map, pdb_folder, lengths, min_residues)
            for pdb_id, pairs in groups
        ]
        for future in futures:
            sizes.update(future.result())

    pp_df["interface_residues"] = pp_df.index.map(sizes)
    pp_df = pp_df[pp_df["interface_residues"] > 0] # remove pairs with no interface residues on the smaller protein's side
    print(f"Computed interface sizes for {len(pp_df)} pairs across {len(groups)} structures", file=log, flush=True)
    return pp_df


if __name__ == "__main__":
    matched_pairs_file = snakemake.input.matched_pairs
    sifts_path = snakemake.input.sifts_file
    pdb_folder = snakemake.input.pdb_folder
    fasta_file = snakemake.input.gene_fasta
    output_file = snakemake.output.interface_annotated
    log_file = snakemake.log[0]

    log = open(log_file, "w")

    pp_df = pd.read_parquet(matched_pairs_file)
    human_proteins = set(pp_df["uniprot_id_bait"]) | set(pp_df["uniprot_id_prey"])
    print(f"Loaded {len(pp_df)} matched pairs over {len(human_proteins)} proteins", file=log, flush=True)
    chain_map, fusions = build_chain_map(sifts_path, human_proteins, log)

    lengths = read_fasta_lengths(fasta_file)
    pp_df = pp_df[pp_df["uniprot_id_bait"].isin(lengths) & pp_df["uniprot_id_prey"].isin(lengths)]
    print(f"{len(pp_df)} pairs with known protein lengths", file=log, flush=True)

    pp_df = set_interface_sizes(pp_df, chain_map, pdb_folder, lengths, log, max_workers=snakemake.threads)
    pp_df.to_parquet(output_file, index=False)
    print("Done.", file=log, flush=True)
    log.close()
