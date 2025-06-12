import re

import pandas
import pandas as pd


def parse_cellosaurus(filename, output_filename):
    rows = []

    with open(filename, "r") as f:
        entry_dict = dict()
        header_passed = False
        for line in f:
            if not header_passed:
                if line[0:2] == "//":
                    header_passed = True
            else:
                try:
                    code, content = line.strip().split("   ")
                    if code == "AC":
                        rows.append(entry_dict)
                        entry_dict = {
                            code: content
                        }
                    elif code == "OX":
                        match = re.search("NCBI_TaxID=([0-9]+)", content)
                        if match:
                            match = match.groups()[0]
                        entry_dict[code] = match

                    elif code == "CA":
                        entry_dict[code] = content
                except  ValueError:
                    continue

    cellosaurus_df = pd.DataFrame(rows)
    cellosaurus_df = cellosaurus_df.rename(
        {
            "AC": "cl_id",
            "OX": "taxon_id",
            "CA": "cl_category"
        },
        axis = 1
    )

    cellosaurus_df.to_csv(output_filename,
        sep="\t",
        index=False
    )


parse_cellosaurus("cellosaurus/cellosaurus.txt",
                  "cellosaurus/cellosaurus_taxonomi_category.csv")
cellosaurus_df = pd.read_csv(
    "cellosaurus/cellosaurus_taxonomi_category.csv",
    sep="\t",
    dtype="str"
)
pubtator_df = pd.read_csv(
    "pubtator/cellline2pubtator3_noblank.csv",
    sep ="\t",
    header = None)
colnames = [
    "pubmed_id",
    "info_type",
    "cl_id",
    "cl_verbose",
    "source"
]
pubtator_df = pubtator_df.rename({i:name for i, name in enumerate(colnames)},axis=1)
pubtator_taxon_df = pubtator_df.merge(cellosaurus_df, on="cl_id")
pub_ids = pubtator_taxon_df[pubtator_taxon_df["taxon_id"] == "9606"]["pubmed_id"].unique()

human_df = pubtator_taxon_df[pubtator_taxon_df["pubmed_id"].isin(pub_ids)]
human_df.to_csv("pubtator/cellline2pubtator3_human.csv", sep="\t", index=None)

pids_cl_count = human_df.groupby(["pubmed_id"])["cl_id"].count()
single_cl_pid_set = set(pids_cl_count[pids_cl_count == 1].index)
single_cl_pid_human_df = human_df[human_df["pubmed_id"].isin(single_cl_pid_set)]
single_cl_pid_human_df.to_csv("pubtator/cellline2pubtator3_human_single_cl.csv", sep="\t", index=None)


