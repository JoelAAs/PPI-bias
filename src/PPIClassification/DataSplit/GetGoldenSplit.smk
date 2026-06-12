import pandas as pd
## https://figshare.com/articles/dataset/PPI_prediction_from_sequence_gold_standard_dataset/21591618/4

rule set_gene_names:
    input:
        gene_names="work_folder/gene_names/gene_names.csv",
        gs_train_pos="data/golden_standard_split/21591618/Intra1_pos_rr.txt",
        gs_train_neg="data/golden_standard_split/21591618/Intra1_neg_rr.txt",
        gs_validation_pos="data/golden_standard_split/21591618/Intra0_pos_rr.txt",
        gs_validation_neg="data/golden_standard_split/21591618/Intra0_neg_rr.txt",
        gs_test_pos="data/golden_standard_split/21591618/Intra2_pos_rr.txt",
        gs_test_neg="data/golden_standard_split/21591618/Intra2_neg_rr.txt"
    output:
        gs_train_pos="work_folder/subsets/train/goldensplit/data_pos.csv",
        gs_train_neg="work_folder/subsets/train/goldensplit/data_neg.csv",
        gs_validation_pos="work_folder/subsets/validation/goldensplit/data_pos.csv",
        gs_validation_neg="work_folder/subsets/validation/goldensplit/data_neg.csv",
        gs_test_pos="work_folder/subsets/test/goldensplit/data_pos.csv",
        gs_test_neg="work_folder/subsets/test/goldensplit/data_neg.csv"
    log:
        "logs/subsets/goldensplit/set_gene_names.log"
    run:
        gene_names = pd.read_csv(input.gene_names,sep="\t")
        for sp_id, gene_id in zip(input[1:], output):
            sp_id_df = pd.read_csv(
                sp_id,sep=" ",header=None)
            sp_id_df.columns = ["uniprot_id_bait", "uniprot_id_prey"]
            sp_id_df = sp_id_df.merge(
                gene_names,left_on="uniprot_id_bait",right_on="uniprot_id")
            del sp_id_df["uniprot_id"]
            sp_id_df = sp_id_df.merge(
                gene_names,left_on="uniprot_id_prey",right_on="uniprot_id",suffixes=("_bait", "_prey"))
            sp_id_df = sp_id_df[["gene_name_bait", "gene_name_prey"]]
            sp_id_df_flipped = sp_id_df.rename(
                {"gene_name_bait":"gene_name_prey", "gene_name_prey":"gene_name_bait"}, axis=1)
            sp_id_df_dir = pd.concat([sp_id_df, sp_id_df_flipped], ignore_index=True) # as they without direction
            sp_id_df_dir.drop_duplicates(inplace=True)

            sp_id_df_dir.to_csv(gene_id, sep="\t", index=False, header=False)

