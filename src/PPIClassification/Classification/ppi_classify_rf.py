import argparse
import datetime

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import randint, uniform
from sklearn.metrics import accuracy_score, classification_report, balanced_accuracy_score
from sklearn.model_selection import ParameterSampler
from sklearn.metrics import f1_score

global RANDOM_STATE


def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length):
    df_pos = pd.read_csv(pos_data_file, sep="\t", usecols=[0, 1], header=None)
    df_negative = pd.read_csv(neg_data_file, sep="\t", usecols=[0, 1], header=None, comment="#")
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
    print("Hyperparameter tuning started")

    param_dist = {
        "n_estimators": randint(400, 2000),
        "max_depth": [None, 8, 10, 12, 16, 20, 24],
        "min_samples_split": randint(10, 500),
        "min_samples_leaf": randint(20, 300),
        "ccp_alpha": [0.0, 1e-5, 1e-4, 1e-3],
        "max_features": [
            "sqrt",
            "log2",
            0.05, 0.1, 0.2, 0.3
        ],
        "bootstrap": [True],
        "max_samples": uniform(0.5, 0.5),
        "min_impurity_decrease": [0.0, 1e-4, 1e-3, 1e-2],
        "class_weight": ["balanced", "balanced_subsample"]
    }

    best_score = -np.inf
    best_model = None
    best_params = None

    for i, params in enumerate(ParameterSampler(param_dist, n_iter=n_iters, random_state=RANDOM_STATE)):
        print(f"{i} of {n_iters} parameter iterations")
        model = RandomForestClassifier(
            **params,
            random_state=RANDOM_STATE,
            n_jobs=n_threads
        )
        s = datetime.datetime.now()
        model.fit(X_train, y_train)
        e = datetime.datetime.now()
        fileout.write("---------------------")
        fileout.write(f"Training took {e - s} for with {params['n_estimators']} estimators on {n_threads} threads\n")
        fileout.write("Current params: " + str(params) + "\n")

        y_test_pred = model.predict(X_validation)
        score = f1_score(y_validation, y_test_pred, average="macro")
        fileout.write(f"Current score: {score} for Validation current n_estimators: {model.n_estimators}\n")

        y_train_pred = model.predict(X_train)
        t_score = f1_score(y_train, y_train_pred, average="macro")
        fileout.write(f"Current score: {t_score} for Train current n_estimators: {model.n_estimators}\n")

        if score > best_score:
            best_score = score
            best_model = model
            best_params = params

    fileout.write("Best validation score: " + str(best_score) + "\n")
    fileout.write("Best params: " + str(best_params) + "\n")

    return best_model, score, params


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
    parser.add_argument("--threads", type=int, default=40, help="")
    parser.add_argument("--randomstate", type=int, default=1234, help="")
    args = parser.parse_args()

    RANDOM_STATE = args.randomstate
    threads = args.threads

    print("Creating embedding dict ... ")
    embed_dict, n_embedding = get_embedding_dict(args.protein_embeddings)

    print("Reading training data ... ")
    X_train, y_train = get_dataset(
        args.train_ppi_data_pos,
        args.train_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    print("Reading validation data ... ")
    X_validate, y_validate = get_dataset(
        args.validation_ppi_data_pos,
        args.validation_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    print("Reading test data ... ")
    X_test, y_test = get_dataset(
        args.test_ppi_data_pos,
        args.test_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    param_file = open(args.params_out, "w")
    _, score, parameters = hyperparameter_tuned_model(X_train, y_train, X_validate, y_validate, threads, param_file, n_iters = 2)

    # DON'T TOUCH UNTIL MIDSOMMAR
    # rfc = RandomForestClassifier(
    #     **parameters,
    #     n_jobs=threads)
    #
    # rfc.fit(
    #     np.vstack((X_train, X_test)),
    #     np.concatenate((y_train, y_test))
    # )
    #
    # y_validate_pred = rfc.predict(X_validate)
    #
    # print("Final test accuracy:", accuracy_score(y_validate, y_validate_pred))
    # print(classification_report(y_test, y_validate_pred))
