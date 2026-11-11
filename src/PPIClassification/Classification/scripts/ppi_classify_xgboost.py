import argparse
from classification_support import hyperparameter_tuned_model
import pandas as pd
from skopt.space import Integer, Real
from xgboost import XGBClassifier
import joblib


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--train_ppi_data_pos", required=True, help="")
    parser.add_argument("--train_ppi_data_neg", required=True, help="")
    parser.add_argument("--protein_embeddings", required=True, help="Path to output csv file")
    parser.add_argument("--params_out", required=True, help="Path to output csv file")
    parser.add_argument("--saved_model", required=True, help="Path to output csv file")
    parser.add_argument("--threads", type=int, default=40, help="")
    parser.add_argument("--randomstate", type=int, default=1234, help="")
    parser.add_argument("--network_type", required=True, type=str, help="Network type, how to setup data")

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
        Integer(3, 10, name="max_depth"),
        Real(0.01, 0.3, prior="log-uniform", name="learning_rate"),
        Integer(1, 50, name="min_child_weight"),
        Real(0.3, 1.0, name="colsample_bytree"),
        Real(0.5, 1.0, name="subsample"),
    ]
    static_parameters = dict(
        n_estimators=100,
        gamma=0,
        use_label_encoder=False,
        eval_metric="logloss",
        verbosity=0,
    )

    param_file = open(args.params_out, "w")
    model, parameters = hyperparameter_tuned_model(
        df_pos, df_negative,
        threads, directed,
        param_file,
        XGBClassifier,
        param_dist,
        static_parameters,
        args.protein_embeddings,
        RANDOM_STATE,
        n_iters=20)

    joblib.dump(model, args.saved_model)
    param_file.close()
