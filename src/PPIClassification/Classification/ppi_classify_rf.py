import argparse
import datetime

import numpy as np
import pandas as pd
from skopt import Optimizer
from skopt.space import Integer, Real, Categorical
from sklearn.ensemble import RandomForestClassifier
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


def hyperparameter_tuned_model(X_train_full, y_train_full, X_validation, y_validation, n_threads, fileout, n_iters=60):
    print("Hyperparameter tuning started", flush=True)

    param_dist = [
        Categorical([None, 8, 10, 12, 16, 20, 24], name="max_depth"),
        Integer(20, 300, name="min_samples_leaf"),
        Categorical(["sqrt", "log2", 0.05, 0.1, 0.2, 0.3], name="max_features"),
        Real(0.1, 0.3, name="max_samples")
    ]

    best_score = np.inf
    best_model = None
    best_params = None
    #best_t = 0.5

    hyper_optimizer = Optimizer(
        dimensions=param_dist,
        base_estimator="RF",
        acq_func="EI",
        random_state=RANDOM_STATE
    )

    for i in range(n_iters):
        s = datetime.datetime.now()
        X_train = X_train_full
        y_train = y_train_full

        params = hyper_optimizer.ask()
        params_dict = dict(zip([d.name for d in param_dist], params))

        model = RandomForestClassifier(
            **params_dict,
            n_estimators=100,
            bootstrap=True,
            random_state=RANDOM_STATE,
            n_jobs=n_threads
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


    param_file = open(args.params_out, "w")
    _, score, parameters = hyperparameter_tuned_model(
        X_train, y_train, X_validate, y_validate, threads, param_file, n_iters=30)

    rfc = RandomForestClassifier(
        **parameters,
        n_estimators=100,
        bootstrap=True,
        random_state=RANDOM_STATE,
        n_jobs=threads)

    rfc.fit(
        np.vstack([X_train,X_validate]),
        np.concatenate([y_train, y_validate])
    )
    joblib.dump(rfc, args.saved_model)
    param_file.close()
    

