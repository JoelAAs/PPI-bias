import pandas as pd

rule count_tested_pairs:
    params:
        id_pattern = config["id_pattern"]
    input:
        pod = f"work_folder{pn}/analysis/POD/POD_{{method}}.csv"
    output:
        pp_counts = f"work_folder{pn}/analysis/POD/summary/POD_{{method}}.csv"
    run:
        pod_df = pd.read_csv(input.pod, sep="\t")
        tested_counts = pod_df.groupby(
            ["n_tested", "n_observed"], as_index=False).size().rename({"size":f"{wildcards.method}_count"},axis=1)

        tested_counts.to_csv(output.pp_counts, sep="\t", index=False)


rule join_counts:
    input:
        counts = expand(f"work_folder{pn}/analysis/POD/summary/POD_{{data}}.csv", data = datasets)
    output:
        all_counts = f"work_folder{pn}/analysis/POD/summary/all.csv"
    run:

        all_dfs = [pd.read_csv(single_count, sep="\t") for single_count in input.counts]
        base_df = all_dfs[0]
        for single_df in all_dfs[1:]:
            base_df = base_df.merge(single_df, by=["n_tested",  "n_observed"], how="full").fillna(0)

