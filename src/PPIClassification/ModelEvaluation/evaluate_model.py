import pandas as pd
import numpy as np
import argparse
import joblib
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


def generate_and_plot_performance(model, X_test, y_test, roc_png):
    y_pred = model.predict_proba(X_test)[:, 1]

    fpr, tpr, _ = roc_curve(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred)
    print("plotting ROC curve", flush=True)
    get_roc_plot((fpr, tpr), roc_auc, roc_png)

    return roc_auc


def get_roc_plot(obs_performance, obs_auc, output_png):
    fpr, tpr = obs_performance

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(fpr, tpr, color='blue', label=f'Observed ROC (AUC={obs_auc:.4f})')
    ax.plot([0, 1], [0, 1], linestyle='--', color='grey', label='Random chance')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curve')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()


def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length, flip_and_double):
    df_pos = pd.read_csv(pos_data_file, sep="\t")[["bait", "prey"]]
    df_negative = pd.read_csv(neg_data_file, sep="\t")[["bait", "prey"]]
    if df_negative.shape[0] > df_pos.shape[0]:
        df_negative = df_negative.sample(df_pos.shape[0], random_state = 1234)
    df_samples = pd.concat([df_pos, df_negative], ignore_index=True)

    if not flip_and_double:
        baits = df_samples.iloc[:, 0].to_numpy()
        prey = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0]
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in baits]
        X[:, embed_length:] = [embedding_dict[p] for p in prey]
        y = np.zeros(n_samples, dtype=np.int8)
        y[:df_pos.shape[0]] = 1

    else:
        protein_a = df_samples.iloc[:, 0].to_numpy()
        protein_b = df_samples.iloc[:, 1].to_numpy()
        n_samples = df_samples.shape[0] * 2
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in np.concatenate([protein_a, protein_b])]
        X[:, embed_length:] = [embedding_dict[p] for p in np.concatenate([protein_b, protein_a])]
        y = np.zeros(n_samples, dtype=np.int8)
        y[:df_pos.shape[0]] = 1
        start_flipped = df_pos.shape[0] + df_negative.shape[0]
        y[start_flipped:start_flipped + df_pos.shape[0]] = 1

    return X, y


def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] - 1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate a trained model on a test set")
    parser.add_argument("--pos_data_file", type=str, required=True)
    parser.add_argument("--neg_data_file", type=str, required=True)
    parser.add_argument("--protein_embeddings_file", type=str, required=True)
    parser.add_argument("--model_file", type=str, required=True)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--plot_roc_png", type=str)
    parser.add_argument("--network_type", required=True, type=str, help="Network type, how to setup data")

    args = parser.parse_args()

    if args.network_type == "undirectional":
        flip_and_double = True
    elif args.network_type == "directional":
        flip_and_double = False
    else:
        raise ValueError(f"{args.network_type} is an invalid network value")


    embedding_dict, embed_length = get_embedding_dict(args.protein_embeddings_file)
    X_test, y_test = get_dataset(args.pos_data_file, args.neg_data_file, embedding_dict, embed_length, flip_and_double)
    model = joblib.load(args.model_file)

    roc_auc = generate_and_plot_performance(model, X_test, y_test, args.plot_roc_png)

    with open(args.output_file, "w") as f:
        f.write(f"ROC: {roc_auc:.4f}\n")
        y_test = y_test.astype(np.int32)
        f.write(f"Samples (pos/neg): {sum(y_test)} / {len(y_test) - sum(y_test)}")
