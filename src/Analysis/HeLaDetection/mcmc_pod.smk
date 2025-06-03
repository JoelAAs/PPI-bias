import pandas as pd

def setup_detection_model():
    

rule estimate_bait_interaction:
    params:
        prior = 1
    input:
        pod_base = "data/absent_0_present_1_selected_N07444_M04547.csv",
        cl_specific_interactions = "data/CL_annotated_bait_prey.csv"
    run:
        PPI_interactions = pd.read_csv(input.cl_specific_interactions, sep="\t")
        baseline_pod = pd.read_csv(input.pod_base).fillna(0)

        PPI_interactions = PPI_interactions[
            (PPI_interactions["gene_name_bait"] != PPI_interactions["gene_name_prey"])
        ]
        PPI_interactions["study_id"] = PPI_interactions[
            ["pubmed_id", "detection_method"]
        ].apply(lambda row: "-".join(map(str, row)), axis=1)
        n_studies = PPI_interactions.groupby("gene_name_bait", as_index=False
        )["study_id"].nunique()
        n_studies = n_studies.rename(
            {
                "study_id": "n_tested"
            }, axis = 1
        )
        n_bait_prey_tests = PPI_interactions.groupby(
            ["gene_name_bait", "gene_name_prey"], as_index=False
        )["study_id"].nunique()
        n_bait_prey_tests = n_bait_prey_tests.rename(
            {
                "study_id": "n_observed"
            }, axis = 1
        )
        n_bait_prey_tests = n_bait_prey_tests[
            n_bait_prey_tests["gene_name_prey"].isin(baseline_pod.columns)
        ]
        all_tested = n_bait_prey_tests.merge(
            n_studies, on="gene_name_bait")