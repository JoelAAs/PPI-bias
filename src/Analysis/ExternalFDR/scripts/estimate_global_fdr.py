import pandas as pd
from estimate_study_fdr import fit_mixture


def main():
    pod_df = pd.read_parquet(
        snakemake.input.pod_file,
        columns=["n_observed", "n_tested"])

    max_iter = snakemake.params.n_em_iterations
    global_estimates = fit_mixture(pod_df, max_iter)

    pd.DataFrame([global_estimates]).to_csv(
        snakemake.output.study_fdr, sep="\t", index=False)


if __name__ == "__main__":
    main()
