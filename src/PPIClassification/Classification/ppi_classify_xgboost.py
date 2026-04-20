import argparse
import datetime

import numpy as np
import pandas as pd
from skopt import Optimizer
from skopt.space import Integer, Real
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, balanced_accuracy_score
from sklearn.metrics import log_loss
import joblib
from sklearn.metrics import precision_recall_curve, auc

global RANDOM_STATE


def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length):
    df_pos = pd.read_csv(pos_data_file, sep="\t")[["bait", "prey"]]
    df_negative = pd.read_csv(neg_data_file, sep="\t")[["bait", "prey"]]
    df_samples = pd.concat([df_pos, df_negative], ignore_index=True)

    baits = df_samples.iloc[:, 0].to_numpy()
    prey = df_samples.iloc[:, 1].to_numpy()

    n_samples = df_samples.shape[0]
    X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
    X[:, :embed_length] = [embedding_dict[b] for b in baits]
    X[:, embed_length:] = [embedding_dict[p] for p in prey]

    y = np.zeros(n_samples, dtype=np.int8)
    y[:df_pos.shape[0]] = 1

    return X, y


def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] - 1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


def hyperparameter_tuned_model(X_train, y_train, X_validation, y_validation, n_threads, fileout, n_iters=10):
    print("Hyperparameter tuning started", flush=True)

    param_dist = [
        Integer(50, 500, prior="log-uniform", name="n_estimators"),
        Integer(3, 10, name="max_depth"),
        Real(0.01, 0.3, prior="log-uniform", name="learning_rate"),
        Integer(1, 50, name="min_child_weight"),
        Real(0.0, 5.0, name="gamma"),
        Real(0.3, 1.0, name="colsample_bytree"),
        Real(0.5, 1.0, name="subsample"),
    ]

    best_score = np.inf
    best_model = None
    best_params = None

    hyper_optimizer = Optimizer(
        dimensions=param_dist,
        base_estimator="RF",
        acq_func="EI",
        random_state=RANDOM_STATE
    )

    for i in range(n_iters):
        s = datetime.datetime.now()

        params = hyper_optimizer.ask()
        params_dict = dict(zip([d.name for d in param_dist], params))

        model = XGBClassifier(
            **params_dict,
            use_label_encoder=False,
            eval_metric="logloss",
            random_state=RANDOM_STATE,
            n_jobs=n_threads,
            verbosity=0
        )
        model.fit(X_train, y_train)
        probs = model.predict_proba(X_validation)[:, 1]
        probs = np.clip(probs, 1e-15, 1 - 1e-15)
        c_log_loss = log_loss(y_validation, probs)

        val_acc = balanced_accuracy_score(y_validation, (probs > 0.5).astype(int))
        hyper_optimizer.tell(params, c_log_loss)

        e = datetime.datetime.now()
        print("------------------------------------------------")
        print(f"{i + 1} iteration of {n_iters} in {e - s}")
        print("Current params: " + str(params))
        print(f"log loss: {c_log_loss}\tAcc: {val_acc} for Validation", flush=True)
        fileout.write("---------------------")
        fileout.write(f"Training took {e - s} using {n_threads} threads\n")
        fileout.write("Current params: " + str(params) + "\n")

        if c_log_loss < best_score:
            best_score = c_log_loss
            best_model = model
            best_params = params

    fileout.write("Best validation score: " + str(best_score) + "\n")
    fileout.write("Best params: " + str(best_params) + "\n")
    best_params = dict(zip([d.name for d in param_dist], best_params))
    return best_model, best_score, best_params


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--train_ppi_data_pos", required=True, help="")
    parser.add_argument("--train_ppi_data_neg", required=True, help="")
    parser.add_argument("--validation_ppi_data_pos", required=True, help="Path to output csv file")
    parser.add_argument("--validation_ppi_data_neg", required=True, help="Path to output csv file")
    parser.add_argument("--test_ppi_data_pos", required=True, help="Path to output csv file")
    parser.add_argument("--test_ppi_data_neg", required=True, help="Path to output csv file")
    parser.add_argument("--protein_embeddings", required=True, help="Path to output csv file")
    parser.add_argument("--params_out", required=True, help="Path to output csv file")
    parser.add_argument("--saved_model", required=True, help="Path to output csv file")
    parser.add_argument("--threads", type=int, default=40, help="")
    parser.add_argument("--randomstate", type=int, default=1234, help="")
    args = parser.parse_args()

    RANDOM_STATE = args.randomstate
    threads = args.threads

    print("Creating embedding dict ... ", flush=True)
    embed_dict, n_embedding = get_embedding_dict(args.protein_embeddings)

    print("Reading training data ... ", flush=True)
    X_train, y_train = get_dataset(
        args.train_ppi_data_pos,
        args.train_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    print("Reading validation data ... ", flush=True)
    X_validate, y_validate = get_dataset(
        args.validation_ppi_data_pos,
        args.validation_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    # print("Reading test data ... ", flush=True)
    # X_test, y_test = get_dataset(
    #     args.test_ppi_data_pos,
    #     args.test_ppi_data_neg,
    #     embed_dict,
    #     n_embedding
    # )

    param_file = open(args.params_out, "w")
    _, score, parameters = hyperparameter_tuned_model(
        X_train, y_train, X_validate, y_validate, threads, param_file, n_iters=60)

    xgb = XGBClassifier(
        **parameters,
        use_label_encoder=False,
        eval_metric="logloss",
        random_state=RANDOM_STATE,
        n_jobs=threads,
        verbosity=0
    )

    xgb.fit(
        np.vstack([X_train, X_validate]),
        np.concatenate([y_train, y_validate])
    )
    joblib.dump(xgb, args.saved_model)

    # probs_test = xgb.predict_proba(X_test)[:, 1]
    # 
    # precision, recall, _ = precision_recall_curve(y_test, probs_test)
    # pr_auc = auc(recall, precision)
    # y_test_pred = (probs_test > 0.5).astype(int)
    # param_file.write("-----------------TEST ACCURACY----------------\n")
    # param_file.write(f"Precision-Recall AUC: {pr_auc:.4f}\n")
    # param_file.write(classification_report(y_test, y_test_pred))
    param_file.close()
