import argparse
import datetime

import numpy as np
import pandas as pd
from skopt import Optimizer
from skopt.space import Integer, Real, Categorical
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, balanced_accuracy_score
from sklearn.metrics import f1_score
import joblib
from sklearn.metrics import precision_recall_curve, auc
from sklearn.dummy import DummyClassifier

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


def hyperparameter_tuned_model(X_train_full, y_train_full, X_validation, y_validation, n_threads, fileout, n_iters=10,
                               max_samples=50000):
    print("Hyperparameter tuning started", flush=True)

    param_dist = [
        Integer(48, 2000, prior="log-uniform", name="n_estimators"),
        Categorical([None, 8, 10, 12, 16, 20, 24], name="max_depth"),
        Integer(10, 500, name="min_samples_split"),
        Integer(20, 300, name="min_samples_leaf"),
        Categorical(["sqrt", "log2", 0.05, 0.1, 0.2, 0.3], name="max_features"),
        Real(0.1, 0.3, name="max_samples")
    ]

    best_score = -np.inf
    best_model = None
    best_params = None
    best_t = 0.5

    hyper_optimizer = Optimizer(
        dimensions=param_dist,
        base_estimator="RF",
        acq_func="EI",
        random_state=RANDOM_STATE
    )

    for i in range(n_iters):
        s = datetime.datetime.now()
        if len(y_train_full) > max_samples:
            rng = np.random.RandomState(RANDOM_STATE + i)
            indices = rng.choice(len(y_train_full), size=max_samples, replace=False)
            X_train = X_train_full[indices]
            y_train = y_train_full[indices]
        else:
            X_train = X_train_full
            y_train = y_train_full

        params = hyper_optimizer.ask()
        params_dict = dict(zip([d.name for d in param_dist], params))

        model = RandomForestClassifier(
            **params_dict,
            bootstrap=True,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=n_threads
        )
        model.fit(X_train, y_train)
        probs = model.predict_proba(X_validation)[:, 1]

        y_val_pred = np.zeros(len(y_validation))
        best_f1 = 0
        current_t = 0.5
        for t in np.linspace(0.1, 0.9, 50):
            preds = (probs > t).astype(int)
            f1 = f1_score(y_validation, preds, average="weighted")
            if f1 > best_f1:
                best_f1 = f1
                current_t = t
                y_val_pred = preds

        val_acc = balanced_accuracy_score(y_validation, y_val_pred)
        hyper_optimizer.tell(params, -best_f1)

        e = datetime.datetime.now()
        print("------------------------------------------------")
        print(f"{i + 1} iteration of {n_iters} in {e - s}")
        print("Current params: " + str(params))
        print(f"F1 score: {best_f1}\t t: {current_t} Acc: {val_acc} for Validation", flush=True)
        fileout.write("---------------------")
        fileout.write(f"Training took {e - s} using {n_threads} threads\n")
        fileout.write("Current params: " + str(params) + "\n")

        if best_f1 > best_score:
            best_score = best_f1
            best_t = current_t
            best_model = model
            best_params = params

    fileout.write("Best validation score: " + str(best_score) + "\n")
    fileout.write("Best params: " + str(best_params) + "\n")
    best_params = dict(zip([d.name for d in param_dist], best_params))
    return best_model, best_t, best_score, best_params



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
    parser.add_argument("--saved_dummy_classifer", required=True, help="Path to output csv file")
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

    print("Reading test data ... ", flush=True)
    X_test, y_test = get_dataset(
        args.test_ppi_data_pos,
        args.test_ppi_data_neg,
        embed_dict,
        n_embedding
    )

    param_file = open(args.params_out, "w")
    _, best_t, score, parameters = hyperparameter_tuned_model(
        X_train, y_train, X_validate, y_validate, threads, param_file, n_iters=20)

    rfc = RandomForestClassifier(
        **parameters,
        bootstrap=True,
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=threads)

    rfc.fit(
        np.vstack((X_train, X_validate)),
        np.concatenate((y_train, y_validate))
    )
    joblib.dump(rfc, args.saved_model)

    probs_test = rfc.predict_proba(X_test)[:, 1]

    precision, recall, _ = precision_recall_curve(y_test, probs_test)
    pr_auc = auc(recall, precision)
    y_test_pred = (probs_test > best_t).astype(int)
    param_file.write("-----------------TEST ACCURACY----------------\n")
    param_file.write(f"Precision-Recall AUC: {pr_auc:.4f}\n")
    param_file.write(f"Selected t: {best_t}\n")
    param_file.write(classification_report(y_test, y_test_pred))
    param_file.close()


    dummy_clf = DummyClassifier(strategy="stratified", random_state=RANDOM_STATE)
    dummy_clf.fit(
        np.vstack((X_train, X_validate)),
        np.concatenate((y_train, y_validate)
        ))
    joblib.dump(dummy_clf, args.saved_dummy_classifer)
    

    