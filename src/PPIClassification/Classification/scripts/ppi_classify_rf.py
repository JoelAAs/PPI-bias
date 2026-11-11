import argparse
import datetime
from classification_support import hyperparameter_tuned_model
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from skopt.space import Integer, Real, Categorical
import joblib

global RANDOM_STATE


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--train_ppi_data_pos", required=True, help="")
    parser.add_argument("--train_ppi_data_neg", required=True, help="")
    parser.add_argument("--protein_embeddings", required=True, help="Path to output csv file")
    parser.add_argument("--params_out", required=True, help="Path to output csv file")
    parser.add_argument("--saved_model", required=True, help="Path to output csv file")
    parser.add_argument("--threads", type=int, default=40, help="")
    parser.add_argument("--randomstate", type=int, default=1234, help="")
    parser.add_argument("--network_type",required=True, type=str, help="Network type, how to setup data")

    args = parser.parse_args()

    RANDOM_STATE = args.randomstate
    threads = args.threads
    if args.network_type == "undirectional":
        directed = False
    elif args.network_type == "directional":
        directed = True
    else:
        raise ValueError(f"{args.network_type} is an invalid network value")
    
    df_pos = pd.read_csv(args.train_ppi_data_pos, sep="\t")[["bait", "prey"]]
    df_negative = pd.read_csv(args.train_ppi_data_neg, sep="\t")[["bait", "prey"]]

    

    param_dist = [
        Categorical([None, 8, 10, 12, 16, 20, 24], name="max_depth"),
        Integer(20, 300, name="min_samples_leaf"),
        Categorical(["sqrt", "log2", 0.05, 0.1, 0.2, 0.3], name="max_features"),
        Real(0.1, 0.3, name="max_samples"),
    ]
    static_parameters = dict(n_estimators=100, bootstrap=True)

    param_file = open(args.params_out, "w")
    model, parameters =  hyperparameter_tuned_model(
        df_pos, df_negative,
        threads, directed,
        param_file,
        RandomForestClassifier,
        param_dist,
        static_parameters,
        args.protein_embeddings,
        RANDOM_STATE,
        n_iters=60)
    
    joblib.dump(model, args.saved_model)
    param_file.close()
    

