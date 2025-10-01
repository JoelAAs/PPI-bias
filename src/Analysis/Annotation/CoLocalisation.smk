import json
import re

import pandas as pd

from localisation_support import *
from mean_distance_support import get_cumulative_sum


def read_localisation_json(filename):
    with open(filename,"r") as f:
        localisation_data = json.load(f)
    return localisation_data


def combine(dict_1, dict_2):
    localisation_keys = set(dict_1) | set(dict_2)
    combined = {k: dict_1.get(k,0) + dict_2.get(k,0) for k in localisation_keys}
    return combined


def get_localisation_per_bait(baits, localisation_files):
    bait_localisations = {b: dict() for b in baits}
    for lf in localisation_files:
        loc_data_dict = read_localisation_json(lf)
        for b in baits:
            bait_localisations[b] = combine(
                bait_localisations[b],
                loc_data_dict.get(b,loc_data_dict["other"]))
        bait_localisations["other"] = combine(
            bait_localisations.get("other", dict()),
            loc_data_dict["other"]
        )

    return bait_localisations


rule method_comparison:
    """
    Compare localisation
    """
    params:
        localisation_csv=config["localisation_file"],
        min_localisation_genes=100
    input:
        multi_method_ms="work_folder/inferred_search_space/aggregated/multi_methods/ms_experimental_wise.csv",
        multi_method_y2h="work_folder/inferred_search_space/aggregated/multi_methods/y2h_experimental_wise.csv"
    output:
        method_diff_localisation="work_folder/inferred_search_space/analysis/localisation/same_localisation_method_diff.csv",
        ms_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_ms.csv",
        y2h_diff_localisation="work_folder/inferred_search_space/analysis/localisation/diff_localisation_y2h.csv"
    run:
        df_ms = pd.read_csv(input.multi_method_ms,sep="\t")
        df_y2h = pd.read_csv(input.multi_method_y2h,sep="\t")

        df_localisation = pd.read_csv(params.localisation_csv,sep="\t")
        localisation_count = df_localisation.groupby("localisation",as_index=False).size()
        keep_localisations = localisation_count[
            localisation_count["size"] > params.min_localisation_genes
            ]["localisation"]
        df_localisation = df_localisation[
            df_localisation["localisation"].isin(keep_localisations)
        ]

        df_ms = add_localisation(df_ms,df_localisation)
        df_ms = df_ms.groupby(['localisation_bait', 'localisation_match'],as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_ms["method"] = "ms"

        df_y2h = add_localisation(df_y2h,df_localisation)
        df_y2h = df_y2h.groupby(['localisation_bait', 'localisation_match'],as_index=False).agg({
            'n_tested': 'sum',
            'n_observed': 'sum'
        })
        df_y2h["method"] = "y2h"

        full_df = pd.concat([df_y2h, df_ms])
        unseen = full_df.groupby("localisation_bait",as_index=False).size()
        localisation_unseen = unseen[unseen["size"] != 4]["localisation_bait"]
        full_df = full_df[~full_df["localisation_bait"].isin(localisation_unseen)]
        full_df["not_observed"] = full_df["n_tested"] - full_df["n_observed"]
        localisations = full_df["localisation_bait"].unique()

        diff_true = fisher_exact_to_df(
            full_df[full_df["localisation_match"] == True],
            localisations,
            "method"
        )
        diff_true.to_csv(output.method_diff_localisation,sep="\t",index=False)

        diff_ms = fisher_exact_to_df(
            full_df[
                full_df["method"] == "ms"
                ].sort_values('localisation_match',ascending=False),
            localisations,
            "localisation_match"
        )
        diff_ms.to_csv(output.ms_diff_localisation,sep="\t",index=False)

        diff_y2h = fisher_exact_to_df(
            full_df[
                full_df["method"] == "y2h"
                ].sort_values('localisation_match',ascending=False),
            localisations,
            "localisation_match"
        )
        diff_y2h.to_csv(output.y2h_diff_localisation,sep="\t",index=False)


rule accumulation_colocalisation:
    params:
        localisation_csv=config["localisation_file"]
    input:
        pod_data="work_folder/analysis/POD/POD_{data}.csv"
    output:
        localisation_annotated="work_folder/analysis/localisation/POD_{data}_localisation.csv",
        localisation_lesser="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_lesser.csv",
        localisation_greater="work_folder/analysis/localisation/cumulative/POD_{data}_localisation_greater.csv"
    run:
        df_localisation = pd.read_csv(params.localisation_csv,sep="\t")
        bait_model = pd.read_csv(input.pod_data,sep="\t")
        proteins_tested = set(bait_model["gene_name_bait"].tolist() + bait_model["gene_name_prey"].tolist())
        df_localisation = df_localisation[df_localisation["gene_name"].isin(proteins_tested)]  # detectable
        n_possible = df_localisation.shape[0]
        random_match = pd.DataFrame([
            [loc, ((df_localisation["localisation"] == loc).sum() - 1) / n_possible] for loc in
            df_localisation["localisation"].unique()
        ],columns=("localisation_prey", "match_probability"))

        bait_model = add_localisation(bait_model,df_localisation)
        bait_model = bait_model.merge(random_match,on="localisation_prey")
        bait_model.to_csv(
            output.localisation_annotated,sep="\t",index=False)

        measurement_columns = ["match_probability", "localisation_match"]
        get_cumulative_sum(
            bait_model,
            value_column="lower_bound_pod",
            cumulative_columns=measurement_columns).to_csv(
            output.localisation_greater,sep="\t",index=False)
        get_cumulative_sum(
            bait_model,
            value_column="upper_bound_pod",
            cumulative_columns=measurement_columns,
            greater=False).to_csv(
            output.localisation_lesser,sep="\t",index=False)


def get_expected_localisations(wc):
    # TODO: HARDCODED and ugly fix, later
    _ = checkpoints.all_methods_filter_out.get(data=wc.data).output[0]
    pod_df = pd.read_csv(f"work_folder/analysis/POD/POD_{wc.data}.csv",sep="\t")
    pids = pod_df["pubmed_id"].unique()
    pids = [p.split(";") for p in pids]
    pids = [item for studies in pids for item in studies]
    pids = set(pids)
    expected = [
        f"work_folder/analysis/localisation/study_match_probability/{pid}.json" for pid in pids
    ]
    return expected


rule get_per_study_localisation:
    params:
        localisation_csv=config["localisation_file"]
    input:
        study="work_folder/inferred_search_space/experimental_method/{pid}_{method}.csv"
    output:
        probability="work_folder/analysis/localisation/study_match_probability/{pid}_{method}.json"
    run:
        study_df = pd.read_csv(input.study,sep="\t")

        df_localisation = pd.read_csv(params.localisation_csv,sep="\t")

        prey_pool = study_df["gene_name_prey"].unique()
        baits = study_df["gene_name_bait"].unique()
        baits = [b for b in baits if b in prey_pool]
        bait_df = df_localisation[df_localisation["gene_name"].isin(baits)]

        prey_localisation = df_localisation[df_localisation["gene_name"].isin(prey_pool)]
        localisation_count = prey_localisation.groupby("localisation",as_index=False).size()
        n_localisations = localisation_count["size"].sum()

        with open(output.probability,"w") as w:
            w.write("{\n")
            loc_tuple = list(localisation_count.itertuples(index=False,name=None))
            for bait in bait_df["gene_name"].unique():
                bait_localisation = bait_df[bait_df["gene_name"] == bait]["localisation"]
                bait_n_localisations = n_localisations - len(bait_localisation)
                if bait_n_localisations == 0:
                    continue
                w.write(f'\t"{bait}": {{\n')
                for i, (localisation, count) in enumerate(loc_tuple):
                    if localisation in bait_localisation:
                        count -= 1
                    w.write(f'\t\t"{localisation}": {str(count / bait_n_localisations)}{"," if i + 1 < len(loc_tuple) else ""}\n')
                w.write("\t},\n")

            w.write(f'\t"other": {{\n')

            for i, (localisation, count) in enumerate(loc_tuple):
                w.write(f'\t\t"{localisation}": {str(count / n_localisations)}{"," if i + 1 < len(loc_tuple) else ""}\n')
            w.write("\t}\n")
            w.write("}\n")

rule get_bait_test_localisation_probability:
    input:
        pod_data="work_folder/analysis/POD/POD_{data}.csv",
        localisations=get_expected_localisations
    output:
        all_probs="work_folder/analysis/localisation/study_match_probability/subsets/{data}_unique_prob.json"
    run:
        unique_probs = pd.read_csv(input.pod_data,sep="\t")[["gene_name_bait", "pubmed_id"]].drop_duplicates()
        unique_probs = unique_probs.groupby("pubmed_id",as_index=False)["gene_name_bait"].unique()
        first = True
        with open(output.all_probs,"w") as w:
            w.write("{\n")
            for _, (pids, baits) in unique_probs.iterrows():
                if not first:
                    w.write(",\n")
                else:
                    first = False
                localisation_files = [
                    f"work_folder/analysis/localisation/study_match_probability/{pid}.json" for pid in pids.split(";")
                ]
                pids_expected = get_localisation_per_bait(baits,localisation_files)
                json_str = json.dumps(pids_expected)
                w.write(f'"{pids}": {json_str}')
            w.write("\n}\n")


def get_study_combination_dict(filename, pubmed_id):
    pattern = '^\"([0-9-_a-zA-Z;]+)\":'
    i = 0
    with open(filename,"r") as f:
        for line in f:
            i += 1
            search = re.search(pattern,line)
            if search:
                current_match = search.groups()[0]
                if pubmed_id == current_match:
                    line = line.strip()
                    if line[-1] == ",":
                        line = line[:-1]
                    prob_dict = json.loads("{" + line + "}")
                    return prob_dict

    raise KeyError(f"{pubmed_id} not in {filename}")

rule annotate_per_study_prob:
    params:
        localisation_csv=config["localisation_file"]
    input:
        pod_data="work_folder/analysis/POD/POD_{data}.csv",
        all_probs="work_folder/analysis/localisation/study_match_probability/subsets/{data}_unique_prob.json"
    output:
        expected_df="work_folder/analysis/localisation/study_match_probability/expected/POD_{data}_expected.csv"
    run:
        localisation_df = pd.read_csv(params.localisation_csv,sep="\t")
        localisation_dict = localisation_df.groupby('gene_name')['localisation'].apply(set).to_dict()

        pod_df = pd.read_csv(input.pod_data,sep="\t")
        pod_df = pod_df.sort_values("pubmed_id")

        pre_pubmed_id = ""
        localisation_prob_dict = dict()
        with open(output.expected_df,"w") as w:
            w.write("\t".join(pod_df.columns.values) + "\tmatch_probability\tlocalisation_match\n")

            for i, row in pod_df.iterrows():

                bait = row["gene_name_bait"]
                prey = row["gene_name_prey"]
                pubmed_id = row["pubmed_id"]
                n_tested = row["n_tested"]
                if pre_pubmed_id != pubmed_id:
                    localisation_prob_dict = get_study_combination_dict(input.all_probs, pubmed_id)
                    pre_pubmed_id = pubmed_id


                match = int(bool(localisation_dict.get(bait,set()) & localisation_dict.get(prey,set())))
                expected = sum(
                    [
                        localisation_prob_dict[pubmed_id].get(
                            bait,localisation_prob_dict[pubmed_id]["other"]
                        ).get(b_localisation,0) for b_localisation in localisation_dict.get(bait,set())
                    ]) / int(n_tested)

                w.write("\t".join(map(str,row.values)) + f"\t{str(expected)}\t{str(match)}\n")


rule accumulation_matched_colocalisation:
    input:
        expected_df="work_folder/analysis/localisation/study_match_probability/expected/POD_{data}_expected.csv"
    output:
        localisation_lesser="work_folder/analysis/localisation/study_match_probability/cumulative/POD_{data}_localisation_lesser.csv",
        localisation_greater="work_folder/analysis/localisation/study_match_probability/cumulative/POD_{data}_localisation_greater.csv"
    run:
        mco_df = pd.read_csv(
            input.expected_df,sep="\t"
        )
        measurement_columns = ["match_probability", "localisation_match"]
        mco_df = mco_df[~((mco_df['match_probability'] == 0) & (mco_df['localisation_match'] == 0))]

        get_cumulative_sum(
            mco_df,
            value_column="lower_bound_pod",
            cumulative_columns=measurement_columns
        ).to_csv(
            output.localisation_greater,
            sep="\t",index=False
        )
        get_cumulative_sum(
            mco_df,
            value_column="upper_bound_pod",
            cumulative_columns=measurement_columns,
            greater=False
        ).to_csv(
            output.localisation_lesser,
            sep="\t",index=False
        )
