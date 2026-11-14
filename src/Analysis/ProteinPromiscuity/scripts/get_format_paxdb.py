
import io
import re
import urllib.request
import zipfile

import pandas as pd




def read_paxdb():
    # pax-db.org wants a named user-agent
    HEADERS = {"User-Agent": "MappingAA"}

    paxdb = pd.read_csv(
        "https://pax-db.org/downloads/6.1/datasets/9606/9606-WHOLE_ORGANISM-integrated.txt",
        sep="\t", comment="#", header=None,
        names=["symbol", "string_external_id", "abundance"],
        storage_options=HEADERS
    )
    
    idmapping = pd.read_csv(
        snakemake.input.idmapping, sep="\t", header=None,
        names=["uniprot_id", "type", "value"]
    )

    string_to_ac = idmapping.loc[idmapping["type"] == "STRING", ["uniprot_id", "value"]] \
        .drop_duplicates("value", keep=False).rename({
            "value": "string_external_id"
        }, axis=1)

    df = paxdb.merge(string_to_ac, on="string_external_id", how="inner")
    df[["uniprot_id", "symbol", "string_external_id", "abundance"]].to_csv(
        snakemake.output.pax_db_data, index=False, sep="\t"
    )


read_paxdb()
