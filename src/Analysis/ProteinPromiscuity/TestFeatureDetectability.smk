import pandas as pd
import numpy as np




rule get_cytocol_vs_localisation:
    input:
        localisation=config["localisation_file"],
    output:
        features="work_folder/analysis/ProteinPromiscuity/features/nucleoplasm-vs-cytocol.csv",
    params:
        reliability=["Enhanced", "Supported", "Approved"],
    run:
        from mygene import MyGeneInfo

        hpa = pd.read_csv(input.localisation, sep="\t")
        hpa = hpa[hpa["Reliability"].isin(params.reliability)]
        main = hpa["Main location"].fillna("")
        anywhere = main + ";" + hpa["Additional location"].fillna("")
        nuc = main.str.contains("Nucleoplasm") & ~anywhere.str.contains("Cytosol")
        cyt = main.str.contains("Cytosol") & ~anywhere.str.contains("Nucleoplasm")
        genes = pd.DataFrame({"gene": hpa["Gene"], "feature": nuc.astype(int)})[
            nuc | cyt
        ]
        hits = MyGeneInfo().querymany(
            genes["gene"].tolist(),
            scopes="ensembl.gene",
            fields="uniprot.Swiss-Prot",
            species="human",
        )
        mapping = (
            pd.DataFrame(
                [
                    (h["query"], h.get("uniprot", {}).get("Swiss-Prot"))
                    for h in hits
                ],
                columns=["gene", "protein_id"],
            )
            .explode("protein_id")
            .dropna()
        )
        features = genes.merge(mapping, on="gene")[
            ["protein_id", "feature"]
        ].drop_duplicates()
        features.drop_duplicates("protein_id", keep=False).to_csv(
            output.features, index=False, sep="\t"
        )


rule test_binary_replication:
    params:
        id_pattern = config["id_pattern"]
    input:
        feature_file="work_folder/analysis/ProteinPromiscuity/features/{feature}.csv",
        pod_file="work_folder/analysis/POD/{network_type}/POD_{dataset}.pq",
    output:
        seperation="work_folder/analysis/ProteinPromiscuity/features/tests/{dataset}_{network_type}_{feature}.tsv",
    log:
        "logs/analysis/ProteinPromiscuity/features/{dataset}_{network_type}_{feature}.log",
    script:
        "scripts/logistic_regression_feature.py"
