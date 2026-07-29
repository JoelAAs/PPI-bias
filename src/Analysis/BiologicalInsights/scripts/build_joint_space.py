"""
Module 4 - joint MS/Y2H pair space. Inner-joins POD_ms and POD_y2h on the unordered
pair, keeping only pairs tested >= concordance_min_tested times in BOTH assays.
"""
import pandas as pd
import pyarrow.dataset as ds


def load_assay(path, bait_col, prey_col, suffix):
    table = ds.dataset(path).to_table(
        columns=[bait_col, prey_col, "n_tested", "n_observed", "lower_bound_pod"]
    )
    df = table.to_pandas()
    n_in = len(df)

    prot_a = df[[bait_col, prey_col]].min(axis=1)
    prot_b = df[[bait_col, prey_col]].max(axis=1)
    df = df.assign(prot_a=prot_a, prot_b=prot_b)
    n_dup = int(df.duplicated(subset=["prot_a", "prot_b"]).sum())

    if n_dup != 0:
        raise ValueError("Undirected network has unmerged edges")

    df = df.rename(columns={
        "n_tested": f"n_tested_{suffix}",
        "n_observed": f"n_observed_{suffix}",
        "lower_bound_pod": f"lower_bound_pod_{suffix}",
    })
    cols = ["prot_a", "prot_b", f"n_tested_{suffix}", f"n_observed_{suffix}",
            f"lower_bound_pod_{suffix}"]
    return df[cols], n_in


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    bait_col = snakemake.params.bait_col
    prey_col = snakemake.params.prey_col
    min_tested = snakemake.params.min_tested

    ms_df, n_ms_in = load_assay(snakemake.input.pod_ms, bait_col, prey_col, "ms")
    y2h_df, n_y2h_in = load_assay(snakemake.input.pod_y2h, bait_col, prey_col, "y2h")
    print(f"MS: {n_ms_in} rows in",
          file=log, flush=True)
    print(f"Y2H: {n_y2h_in} rows in",
          file=log, flush=True)

    joint = ms_df.merge(y2h_df, on=["prot_a", "prot_b"], how="inner")
    n_joint_all = len(joint)
    print(f"{n_joint_all} pairs tested in both assays (any n_tested)", file=log, flush=True)

    joint = joint[
        (joint["n_tested_ms"] >= min_tested) & (joint["n_tested_y2h"] >= min_tested)
    ].reset_index(drop=True)
    print(f"{len(joint)}/{n_joint_all} pairs kept after requiring n_tested >= "
          f"{min_tested} in both assays", file=log, flush=True)

    joint.to_parquet(snakemake.output.joint_space, index=False)
    print("Done.", file=log, flush=True)
    log.close()
