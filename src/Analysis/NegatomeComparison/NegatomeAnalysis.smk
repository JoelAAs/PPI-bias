import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds


rule negatome_comparison:
    """
    Compare if Negatome2.0 is comparable to high confidence non-interactors
    """
    input:
        pod_data="work_folder/analysis/POD/undirectional/POD_flat.pq",
        negatome2="data/negatome2.tsv"
    output:
        table="work_folder/analysis/neg2compare/negatome2.txt"
    log:
        "logs/analysis/neg2compare/negatome2.log"
    run:
        neg2_df = pd.read_csv(input.negatome2, sep="\t")
        protein_a = neg2_df["ProteinA"].astype(str)
        protein_b = neg2_df["ProteinB"].astype(str)
        # Undirectional: order the two accessions so A-B and B-A collapse to one key.
        neg2_pairs = {
            f"{a}_{b}" if a < b else f"{b}_{a}" for a, b in zip(protein_a, protein_b)
        }
        # POD_flat.pq is ~100M rows. A row can only match a negatome pair if both of its
        # proteins are negatome proteins, so push that down to the parquet reader: it prunes
        # row groups on column statistics and never materialises the rest. The pair key is
        # then built on the few thousand survivors instead of on every row.
        neg2_proteins = pa.array(
            sorted(set(protein_a) | set(protein_b)), type=pa.string()
        )
        candidates = ds.dataset(input.pod_data, format="parquet").to_table(
            columns=["uniprot_id_bait", "uniprot_id_prey", "n_observed"],
            filter=ds.field("uniprot_id_bait").isin(neg2_proteins)
            & ds.field("uniprot_id_prey").isin(neg2_proteins),
        ).to_pandas()

        bait = candidates["uniprot_id_bait"].astype(str)
        prey = candidates["uniprot_id_prey"].astype(str)
        bait_first = bait < prey
        pair_id = bait.where(bait_first, prey) + "_" + prey.where(bait_first, bait)

        joined = candidates[pair_id.isin(neg2_pairs)]
        n_observed = int((joined["n_observed"] != 0).sum())
        n_not_observed = int((joined["n_observed"] == 0).sum())

        with open(output.table, "w") as w:
            print(f"negatome2 rows:\t{len(neg2_df)}\n")
            print(f"negatome2 unique undirectional pairs:\t{len(neg2_pairs)}\n")
            print(f"joined rows, n_observed != 0:\t{n_observed}\n")
            print(f"joined rows, n_observed == 0:\t{n_not_observed}\n")
            print(f"joined rows, total:\t{n_observed + n_not_observed}\n")
