"""
Module 2 - disorder content (BiologicalInsights plan §7.4): fraction of residues in
disordered regions, restricted to hub + reference proteins only (not the whole
proteome - this is a per-accession API call, kept small deliberately).

Prefer MobiDB-lite consensus (content_fraction from the precomputed MobiDB API
record); proteins absent from MobiDB fall back to UniProt's own ft_region entries
whose description contains "Disordered" (fraction of residues covered by such
regions). Source is recorded per protein - never mixed silently.
"""
import json
import re
import time

import pandas as pd
import requests

BATCH_SIZE = 100
BATCH_DELAY_S = 0.2
MAX_RETRIES = 4
REGION_RE = re.compile(r'REGION (\d+)\.\.(\d+); /note="Disordered"')


def _get_with_retries(url, params, log=None):
    """MobiDB resets the connection under rapid back-to-back batch requests
    (observed live) - retry with exponential backoff rather than failing the whole
    rule on one transient batch."""
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.ok:
                return resp
        except requests.exceptions.RequestException as exc:
            if log is not None:
                print(f"  batch request failed (attempt {attempt + 1}/{MAX_RETRIES}): "
                      f"{exc!r}", file=log, flush=True)
        time.sleep(2 ** attempt)
    return None


def fetch_mobidb_batch(accessions, log=None):
    """MobiDB's batch endpoint returns newline-delimited JSON (one record per
    accession), not a JSON array - verified against a live query."""
    out = {}
    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i:i + BATCH_SIZE]
        resp = _get_with_retries(
            "https://mobidb.org/api/download",
            {"acc": ",".join(batch), "format": "json"},
            log=log,
        )
        time.sleep(BATCH_DELAY_S)
        if resp is None:
            if log is not None:
                print(f"  giving up on batch {i // BATCH_SIZE} ({len(batch)} accessions) "
                      f"after {MAX_RETRIES} attempts - left missing, not imputed",
                      file=log, flush=True)
            continue
        for line in resp.text.strip().split("\n"):
            if not line:
                continue
            rec = json.loads(line)
            disorder = rec.get("prediction-disorder-mobidb_lite")
            if disorder and "content_fraction" in disorder:
                out[rec["acc"]] = disorder["content_fraction"]
    return out


def fallback_uniprot_region_fraction(accessions, uniprot_df):
    out = {}
    sub = uniprot_df.set_index("uniprot_id").reindex(accessions)
    for acc, row in sub.iterrows():
        length = row.get("length")
        region_text = row.get("ft_region")
        if pd.isna(length) or length <= 0 or not isinstance(region_text, str):
            continue
        covered = 0
        for start, end in REGION_RE.findall(region_text):
            covered += int(end) - int(start) + 1
        if covered > 0:
            out[acc] = min(covered / length, 1.0)
    return out


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")

    hub_df = pd.read_csv(snakemake.input.hub_set, sep="\t")
    reference_df = pd.read_csv(snakemake.input.reference_set, sep="\t")
    uniprot_df = pd.read_csv(snakemake.input.uniprot_annotation, sep="\t")

    accessions = sorted(set(hub_df["uniprot_id"]) | set(reference_df["uniprot_id"]))
    print(f"Fetching disorder content for {len(accessions)} hub+reference proteins",
          file=log, flush=True)

    mobidb_vals = fetch_mobidb_batch(accessions, log=log)
    print(f"{len(mobidb_vals)}/{len(accessions)} found in MobiDB-lite (primary source)",
          file=log, flush=True)

    missing = [a for a in accessions if a not in mobidb_vals] 
    fallback_vals = fallback_uniprot_region_fraction(missing, uniprot_df)
    print(f"{len(fallback_vals)}/{len(missing)} of the MobiDB-missing proteins recovered "
          f"from UniProt ft_region 'Disordered' entries (fallback source)",
          file=log, flush=True)

    n_no_data = len(accessions) - len(mobidb_vals) - len(fallback_vals)
    print(f"{n_no_data}/{len(accessions)} proteins have no disorder data from either "
          f"source - left missing, not imputed", file=log, flush=True)

    rows = []
    for acc in accessions:
        if acc in mobidb_vals:
            rows.append({"uniprot_id": acc, "disorder_fraction": mobidb_vals[acc],
                          "disorder_source": "mobidb_lite"})
        elif acc in fallback_vals:
            rows.append({"uniprot_id": acc, "disorder_fraction": fallback_vals[acc],
                          "disorder_source": "uniprot_ft_region"})
        else:
            rows.append({"uniprot_id": acc, "disorder_fraction": None,
                          "disorder_source": "none"})

    # Explicit columns so an empty `rows` (hub set was empty upstream) still writes a
    # header - pd.DataFrame([]) has zero columns, which makes the downstream
    # pd.read_csv() fail with EmptyDataError instead of hitting that script's own
    # empty-hub-set handling.
    columns = ["uniprot_id", "disorder_fraction", "disorder_source"]
    pd.DataFrame(rows, columns=columns).to_csv(snakemake.output.disorder, sep="\t", index=False)
    print("Done.", file=log, flush=True)
    log.close()
