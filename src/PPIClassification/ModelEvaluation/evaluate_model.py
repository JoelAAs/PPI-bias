import pandas as pd
import numpy as np
import sys
import joblib
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


def generate_and_plot_performance(model, X_test, y_test, roc_png):
    y_pred = model.predict_proba(X_test)[:, 1]

    fpr, tpr, _ = roc_curve(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred)
    print("plotting ROC curve", flush=True)
    get_roc_plot((fpr, tpr), roc_auc, roc_png)

    y_pred_label = (y_pred >= 0.5).astype(np.int8)
    pos_mask = y_test == 1
    neg_mask = y_test == 0
    pos_accuracy = (y_pred_label[pos_mask] == y_test[pos_mask]).mean()
    neg_accuracy = (y_pred_label[neg_mask] == y_test[neg_mask]).mean()

    return roc_auc, pos_accuracy, neg_accuracy, y_pred


def get_roc_plot(obs_performance, obs_auc, output_png):
    fpr, tpr = obs_performance

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(fpr, tpr, color="blue", label=f"Observed ROC (AUC={obs_auc:.4f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", label="Random chance")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()


def get_dataset(
    pos_data_file, neg_data_file, embedding_dict, embed_length, flip_and_double
):
    df_pos = pd.read_csv(pos_data_file, sep="\t")[["bait", "prey"]]
    df_negative = pd.read_csv(neg_data_file, sep="\t")[["bait", "prey"]]
    df_samples = pd.concat([df_pos, df_negative], ignore_index=True)

    if not flip_and_double:
        baits = df_samples.iloc[:, 0].to_numpy()
        prey = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0]
        pairs = df_samples
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in baits]
        X[:, embed_length:] = [embedding_dict[p] for p in prey]
        y = np.zeros(n_samples, dtype=np.int8)
        y[: df_pos.shape[0]] = 1

    else:
        protein_a = df_samples.iloc[:, 0].to_numpy()
        protein_b = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0] * 2
        flipped = df_samples.rename(columns={"bait": "prey", "prey": "bait"})
        pairs = pd.concat([df_samples, flipped], ignore_index=True)
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

    return X, y, pairs


def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] - 1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


def centred_unit_embeddings(embedding_dict):
    # TODO: Check
    # as in within_pair_similarity_auc.py: centre over all proteins, unit length -> cosine = dot product
    ids = list(embedding_dict)
    E = np.vstack([embedding_dict[i] for i in ids]).astype(np.float64)
    E -= E.mean(axis=0)
    E /= np.linalg.norm(E, axis=1, keepdims=True)
    return dict(zip(ids, E))


def adjusted_auc(y, y_pred, pairs, unit_emb, n_bins=5):
    # TODO: Check
    sim = np.array([unit_emb[b] @ unit_emb[p] for b, p in zip(pairs["bait"], pairs["prey"])])
    bins = pd.qcut(sim, n_bins, labels=False)
    aucs = [roc_auc_score(y[bins == b], y_pred[bins == b]) for b in range(n_bins)]
    return np.average(aucs, weights=[(y[bins == b] == 1).sum() for b in range(n_bins)])


if __name__ == "__main__":
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")

    network_type = snakemake.wildcards.network_type
    if network_type == "undirectional":
        flip_and_double = True
    elif network_type == "directional":
        flip_and_double = False
    else:
        raise ValueError(f"{network_type} is an invalid network value")

    embedding_dict, embed_length = get_embedding_dict(
        snakemake.input.protein_embeddings
    )
    X_test, y_test, pairs_test = get_dataset(
        snakemake.input.test_pos,
        snakemake.input.test_neg,
        embedding_dict,
        embed_length,
        flip_and_double,
    )
    X_test_no, y_test_no, pairs_test_no = get_dataset(
        snakemake.input.test_pos,
        snakemake.input.test_neg_no,
        embedding_dict,
        embed_length,
        flip_and_double,
    )
    model = joblib.load(snakemake.input.saved_model)

    roc_auc, pos_accuracy, neg_accuracy, y_pred = generate_and_plot_performance(
        model, X_test, y_test, snakemake.output.roc_png
    )
    roc_auc_no, _, neg_accuracy_no, y_pred_no = generate_and_plot_performance(
        model, X_test_no, y_test_no, snakemake.output.roc_png_no
    )
    unit_emb = centred_unit_embeddings(embedding_dict)
    roc_auc_adj = adjusted_auc(y_test, y_pred, pairs_test, unit_emb)
    roc_auc_adj_no = adjusted_auc(y_test_no, y_pred_no, pairs_test_no, unit_emb)
    pairs_test.assign(label=y_test, score=y_pred).to_csv(snakemake.output.pred, sep="\t", index=False)
    pairs_test_no.assign(label=y_test_no, score=y_pred_no).to_csv(snakemake.output.pred_no, sep="\t", index=False)
    n_samples_file = lambda path: pd.read_csv(path, sep="\t", usecols=["bait"]).shape[0]

    samples_train_pos = n_samples_file(snakemake.input.train_pos)
    samples_train_neg = n_samples_file(snakemake.input.train_neg)
    

    with open(snakemake.output.metrics, "w") as f:
        f.write(f"ROC: {roc_auc:.4f}\n")
        f.write(f"Accuracy interaction: {pos_accuracy:.4f}\n")
        f.write(f"Accuracy non-interaction: {neg_accuracy:.4f}\n")
        y_test = y_test.astype(np.int32)
        f.write(f"Samples (pos/neg): {sum(y_test)} / {len(y_test) - sum(y_test)}\n")
        f.write(f"ROC NO: {roc_auc_no:.4f}\n")
        f.write(f"Accuracy non-interaction NO: {neg_accuracy_no:.4f}\n")
        f.write(f"Train Samples (pos/neg): {samples_train_pos}/{samples_train_neg}\n")
        y_test_no = y_test_no.astype(np.int32)
        f.write(
            f"Samples NO (pos/neg): {sum(y_test_no)} / {len(y_test_no) - sum(y_test_no)}\n"
        )
        f.write(f"ROC adjusted: {roc_auc_adj:.4f}\n")
        f.write(f"ROC adjusted NO: {roc_auc_adj_no:.4f}\n")
