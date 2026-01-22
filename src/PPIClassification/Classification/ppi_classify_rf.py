import argparse
import datetime

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import randint
from sklearn.metrics import accuracy_score, classification_report, balanced_accuracy_score
from sklearn.model_selection import ParameterSampler

global RANDOM_STATE


def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length):
    df_pos = pd.read_csv(pos_data_file, sep="\t", usecols=[0, 1])
    df_negative = pd.read_csv(neg_data_file, sep="\t", usecols=[0, 1])
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
        "n_estimators": randint(200, 800),
        "max_depth": randint(4, 12),
        "min_samples_split": randint(20, 200),
        "min_samples_leaf": randint(10, 100),
        "max_features": "sqrt",
        "class_weight": "balanced"
    }

    best_score = -np.inf
    best_model = None
    best_params = None

    for params in ParameterSampler(param_dist, n_iter=n_iters, random_state=RANDOM_STATE):

        model = RandomForestClassifier(
            **params,
            random_state=RANDOM_STATE,
            n_jobs=n_threads
        )
        s = datetime.datetime.now()
        model.fit(X_train, y_train)
        e = datetime.datetime.now()
        fileout.write(f"Training took {e-s} for with {params['n_estimators']} estimators on {n_threads} threads")

        y_test_pred = model.predict(X_validation)
        score = balanced_accuracy_score(y_validation, y_test_pred)
        fileout.write(f"Current score: {score} for Validation current n_estimators: {model.n_estimators}")

        y_train_pred = model.predict(X_train)
        t_score = balanced_accuracy_score(y_train, y_train_pred)
        fileout.write(f"Current score: {t_score} for Train current n_estimators: {model.n_estimators}")

        if score > best_score:
            best_score = score
            best_model = model
            best_params = params

    fileout.write("Best validation score:", best_score)
    fileout.write("Best params:", best_params)

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

    RANDOM_STATE=args.randomstate
    threads = args.threads

    print("Creating embedding dict ... ")
    embed_dict, n_embedding =  get_embedding_dict(args.protein_embeddings)

    print("Reading training data ... ")
    X_train, y_train = get_dataset(
        args.train_ppi_data_pos,
        args.train_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    print("Reading training data ... ")
    X_validate, y_validate = get_dataset(
        args.validate_ppi_data_pos,
        args.validate_ppi_data_neg,
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
    _, score, parameters = hyperparameter_tuned_model(X_train, y_train, X_validate, y_validate, threads, param_file)

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
