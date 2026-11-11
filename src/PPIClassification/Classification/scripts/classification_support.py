import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skopt import Optimizer
from sklearn.metrics import balanced_accuracy_score, log_loss

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "DataSplit" / "scripts"))
from sample_balance_multi_network_functions import degree_balace_edges


def get_dataset(df_pos, df_negative, embedding_dict, embed_length, flip_and_double):
    df_samples = pd.concat([df_pos, df_negative], ignore_index=True)

    if not flip_and_double:
        baits = df_samples.iloc[:, 0].to_numpy()
        prey = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0]
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in baits]
        X[:, embed_length:] = [embedding_dict[p] for p in prey]
        y = np.zeros(n_samples, dtype=np.int8)
        y[: df_pos.shape[0]] = 1

    else:
        protein_a = df_samples.iloc[:, 0].to_numpy()
        protein_b = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0] * 2
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [
            embedding_dict[b] for b in np.concatenate([protein_a, protein_b])
        ]
        X[:, embed_length:] = [
            embedding_dict[p] for p in np.concatenate([protein_b, protein_a])
        ]
        y = np.zeros(n_samples, dtype=np.int8)
        y[: df_pos.shape[0]] = 1
        start_flipped = df_pos.shape[0] + df_negative.shape[0]
        y[start_flipped : start_flipped + df_pos.shape[0]] = 1

    return X, y


def validation_train_split(df_train_pos, df_train_neg, directed, random_state, split=0.2):
    rng = np.random.default_rng(random_state)
    if directed:
        bait_proteins = df_train_pos["bait"].unique()
        prey_proteins = df_train_pos["prey"].unique()
        val_bait = rng.choice(
            bait_proteins, size=int(split * len(bait_proteins)), replace=False
        )
        val_prey = rng.choice(
            prey_proteins, size=int(split * len(prey_proteins)), replace=False
        )
        in_val = lambda df: df["bait"].isin(val_bait) & df["prey"].isin(val_prey)
        in_train = lambda df: ~df["bait"].isin(val_bait) & ~df["prey"].isin(val_prey)
    else:
        proteins = np.union1d(df_train_pos["bait"], df_train_pos["prey"])
        val = rng.choice(proteins, size=int(split * len(proteins)), replace=False)
        in_val = lambda df: df["bait"].isin(val) & df["prey"].isin(val)
        in_train = lambda df: ~df["bait"].isin(val) & ~df["prey"].isin(val)

    selected_train_pos, selected_train_neg, _ = degree_balace_edges(
        df_train_pos[in_train(df_train_pos)],
        df_train_neg[in_train(df_train_neg)],
        0.95,
        directed,
        seed=random_state,
    )
    selected_val_pos, selected_val_neg, _ = degree_balace_edges(
        df_train_pos[in_val(df_train_pos)],
        df_train_neg[in_val(df_train_neg)],
        0.95,
        directed,
        seed=random_state,
    )
    return selected_train_pos, selected_train_neg, selected_val_pos, selected_val_neg


def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] - 1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


def hyperparameter_tuned_model(
    df_train_pos,
    df_train_neg,
    n_threads,
    directed,
    fileout,
    model_function,
    param_dist,
    static_parameters,
    protein_embeddings_file,
    random_state,
    n_iters=60,
):
    print("Hyperparameter tuning started", flush=True)

    best_score = np.inf
    best_model = None
    best_params = None

    hyper_optimizer = Optimizer(
        dimensions=param_dist,
        base_estimator="RF",
        acq_func="EI",
        random_state=random_state,
    )

    embed_dict, n_embedding = get_embedding_dict(protein_embeddings_file)

    print("Reading training data ... ", flush=True)
    n_full = df_train_pos.shape[0] + df_train_neg.shape[0]
    for i in range(n_iters):
        s = datetime.datetime.now()
        pos_train_edges, neg_train_edges, pos_val_edges, neg_val_edges = (
            validation_train_split(df_train_pos, df_train_neg, directed, random_state + i)
        )
        X_train, y_train = get_dataset(
            pos_train_edges,
            neg_train_edges,
            embed_dict,
            n_embedding,
            flip_and_double=(not directed),
        )
        X_validation, y_validation = get_dataset(
            pos_val_edges,
            neg_val_edges,
            embed_dict,
            n_embedding,
            flip_and_double=(not directed),
        )

        
        params = hyper_optimizer.ask()
        params_dict = dict(zip([d.name for d in param_dist], params))
        param_dict = {**params_dict, **static_parameters}

        model = model_function(
            **param_dict,
            random_state=random_state,
            n_jobs=n_threads,
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
        n_train = pos_train_edges.shape[0] + neg_train_edges.shape[0]
        n_val = pos_val_edges.shape[0] + neg_val_edges.shape[0]
        fileout.write("---------------------\n")
        fileout.write(
            f"train: {n_train} ({100 * n_train / n_full:.1f}%), "
            f"validation: {n_val} ({100 * n_val / n_full:.1f}%) "
            f"of {n_full} edges\n"
        )
        fileout.write(f"Training took {e - s} using {n_threads} threads\n")
        fileout.write("Current params: " + str(params) + "\n")

        if c_log_loss < best_score:
            best_score = c_log_loss
            best_model = model
            best_params = param_dict

    fileout.write("Best validation score: " + str(best_score) + "\n")
    fileout.write("Best params: " + str(best_params) + "\n")

    model = model_function(**best_params, random_state=random_state, n_jobs=n_threads)

    X_train_full, y_train_full = get_dataset(
        df_train_pos,
        df_train_neg,
        embed_dict,
        n_embedding,
        flip_and_double=(not directed),
    )

    model.fit(X_train_full, y_train_full)

    return model, best_params
