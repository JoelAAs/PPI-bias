from scipy.stats import beta
import datetime


def pod_input(wc):
    if wc.data == "cell_line":
        return "work_folder/inferred_search_space/aggregated/cell_line/cell_line_experimental_wise.csv"
    return f"work_folder/inferred_search_space/aggregated/methods/{wc.data}_experimental_wise.csv"


checkpoint all_methods_filter_out:
    """
    Get upper and lower bound probability of detection given test/observations of each bait-prey combination.
    For data == "cell_line" pairs are kept per CVCL and the prior is estimated per CVCL; other datasets
    have CVCL == NaN for every row, so this reduces to a single global prior.
    """
    params:
        pseudo_n=config["pseudo_n"],
        id_pattern=config["id_pattern"],
        require_fasta=config["require_fasta"]
    input:
        method_aggregate=pod_input,
        fasta="work_folder/protein_sequences/uniprot_canonical.fasta"
    output:
        full_detection="work_folder/analysis/POD/{network_type}/POD_{data}.pq"
    log:
        "logs/analysis/POD/{network_type}/POD_{data}.log"
    run:
        bait_col = f"{params.id_pattern}_bait"
        prey_col = f"{params.id_pattern}_prey"

        fasta_ids = set()
        with open(input.fasta) as _fasta:
            for line in _fasta:
                if line.startswith(">"):
                    parts = line.split("|")
                    if len(parts) >= 2:
                        fasta_ids.add(parts[1])

        inferred_negative_df = pd.read_csv(
            input.method_aggregate,
            sep="\t",
            dtype={
                bait_col: "string",
                prey_col: "string"
            }
        )

        if params.require_fasta:
            inferred_negative_df = inferred_negative_df[
                inferred_negative_df[bait_col].isin(fasta_ids) &
                inferred_negative_df[prey_col].isin(fasta_ids)
            ].copy()

        inferred_negative_df = inferred_negative_df[
            inferred_negative_df[bait_col] != inferred_negative_df[prey_col]
            ].copy()

        # grouping key for cell line, "" for non-cell-line datasets (CVCL is NaN there)
        inferred_negative_df["_cl"] = inferred_negative_df["CVCL"].fillna("").astype(str)

        if wildcards.network_type == "undirectional":
            s = datetime.datetime.now()
            n_before = len(inferred_negative_df)
            bait_arr = inferred_negative_df[bait_col].to_numpy(dtype=object)
            prey_arr = inferred_negative_df[prey_col].to_numpy(dtype=object)
            inferred_negative_df[bait_col] = np.where(bait_arr <= prey_arr, bait_arr, prey_arr)
            inferred_negative_df[prey_col] = np.where(bait_arr <= prey_arr, prey_arr, bait_arr)
            inferred_negative_df["id_var"] = inferred_negative_df[bait_col] + "_" + inferred_negative_df[prey_col]
            orig_cols = inferred_negative_df.columns.tolist()

            key = ["id_var", "_cl"]
            to_merge = inferred_negative_df.duplicated(key, keep=False)  # only pairs seen both ways need aggregating
            merged = (
                inferred_negative_df[to_merge]
                .groupby(key, sort=False)
                .agg({
                    "n_tested": "sum",
                    "n_observed": "sum",
                    "pubmed_id": lambda x: ";".join({p for ids in x for p in ids.split(";")}),
                })
            )
            # keep first-occurrence row order so pair_id matches a full groupby(sort=False)
            inferred_negative_df = inferred_negative_df[~inferred_negative_df.duplicated(key, keep="first")]
            inferred_negative_df = inferred_negative_df.set_index(key)
            inferred_negative_df.update(merged)
            inferred_negative_df = inferred_negative_df.reset_index()[orig_cols]
            print(f"joined {n_before - len(inferred_negative_df)} out of {n_before} rows in {(datetime.datetime.now() - s).total_seconds()} seconds", flush=True)

        inferred_negative_df["n_observed"] = inferred_negative_df["n_observed"].astype(int)
        inferred_negative_df["n_tested"] = inferred_negative_df["n_tested"].astype(int)

        # P(interaction | tested), per cell line (a single global group for non-cell-line data)
        pod_by_cl = inferred_negative_df.groupby("_cl")[["n_observed", "n_tested"]].sum()
        pod_by_cl = pod_by_cl["n_observed"] / pod_by_cl["n_tested"] # if not cell line it will have one value so

        alpha_prior = params.pseudo_n * inferred_negative_df["_cl"].map(pod_by_cl)
        beta_prior = params.pseudo_n - alpha_prior
        del inferred_negative_df["_cl"]

        inferred_negative_df["alpha_post"] = alpha_prior + inferred_negative_df["n_observed"]
        inferred_negative_df["beta_post"] = beta_prior + inferred_negative_df["n_tested"] - inferred_negative_df[
            "n_observed"]

        inferred_negative_df["p"] = inferred_negative_df["alpha_post"] / (
                inferred_negative_df["alpha_post"] + inferred_negative_df["beta_post"])

        inferred_negative_df["lower_bound_pod"] = beta.ppf(0.025,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["upper_bound_pod"] = beta.ppf(0.975,
            inferred_negative_df["alpha_post"],inferred_negative_df["beta_post"])
        inferred_negative_df["pair_id"] = range(inferred_negative_df.shape[0])

        inferred_negative_df.to_parquet(
            output.full_detection,
            index=False
        )
