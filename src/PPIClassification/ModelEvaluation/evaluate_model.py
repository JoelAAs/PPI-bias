import pandas as pd 
import numpy as np 
import argparse
import joblib
from sklearn.metrics import precision_recall_curve, auc, roc_curve
from sklearn.dummy import DummyClassifier


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

def get_dummy_pr_distribution(y_test, pred_dummy, n=1000):
    for _ in range(n):
        pred_dummy = np.random.permutation(pred_dummy)
        precision, recall, _ = precision_recall_curve(y_test, pred_dummy)
        pr_auc_dummy = auc(recall, precision)
        yield pr_auc_dummy


def evaluate_model(model, dummy_model, X_test, y_test):
    pr_auc, pr_auc_dummy, roc_auc, roc_auc_dummy = None, None, None, None
    if hasattr(model, "predict_proba"):
        probs_test = model.predict_proba(X_test)[:, 1]
        precision, recall, _ = precision_recall_curve(y_test, probs_test)
        probs_test_dummy = dummy_model.predict_proba(X_test)[:, 1]

        # PR AUC
        pr_auc = auc(recall, precision) 
        pr_auc_dummy_dist = get_dummy_pr_distribution(y_test, probs_test_dummy)
        pr_auc_dummy = np.mean(list(pr_auc_dummy_dist))

        # PR AUC negative class
        y_test_neg = 1 - y_test
        probs_test_neg = 1 - probs_test
        probs_test_dummy_neg = 1 - probs_test_dummy

        precision_neg, recall_neg, _ = precision_recall_curve(
            y_test_neg, probs_test_neg
        )
        pr_auc_neg = auc(recall_neg, precision_neg)

        pr_auc_dummy_neg_dist = get_dummy_pr_distribution(y_test_neg, probs_test_dummy_neg)
        pr_auc_dummy_neg = np.mean(list(pr_auc_dummy_neg_dist))
        
        # ROC AUC
        fpr, tpr, thresholds = roc_curve(y_test, probs_test)
        roc_auc = auc(fpr, tpr)
        fpr_dummy, tpr_dummy, _ = roc_curve(y_test, probs_test_dummy)
        roc_auc_dummy = auc(fpr_dummy, tpr_dummy)
        
    return pr_auc, pr_auc_dummy, pr_auc_neg, pr_auc_dummy_neg, roc_auc, roc_auc_dummy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate a trained model on a test set")
    parser.add_argument("--pos_data_file", type=str, required=True, help="Path to the positive samples file")
    parser.add_argument("--neg_data_file", type=str, required=True, help="Path to the negative samples file")
    parser.add_argument("--protein_embeddings_file", type=str, required=True, help="Path to the protein embeddings file")
    parser.add_argument("--model_file", type=str, required=True, help="Path to the saved model file")
    parser.add_argument("--dummy_baseline_file", type=str, required=True, help="Path to save dummy baseline results")
    parser.add_argument("--output_file", type=str, default="evaluation_results.txt", help="File to save evaluation results")   
    args = parser.parse_args()

    embedding_dict, embed_length = get_embedding_dict(args.protein_embeddings_file)
    X_test, y_test = get_dataset(args.pos_data_file, args.neg_data_file, embedding_dict, embed_length)

    model = joblib.load(args.model_file)
    dummy_model = joblib.load(args.dummy_baseline_file) 
    
    with open(args.output_file, "w") as f:
        pr_auc, pr_auc_dummy, pr_auc_neg, pr_auc_dummy_neg, roc_auc, roc_auc_dummy = evaluate_model(model, dummy_model, X_test, y_test)
        f.write(f"PR AUC: {pr_auc:.4f}\n")
        f.write(f"PR AUC (Dummy): {pr_auc_dummy:.4f}\n")
        f.write(f"PR NEG AUC: {pr_auc_neg:.4f}\n")
        f.write(f"PR NEG AUC (Dummy): {pr_auc_dummy_neg:.4f}\n")
        f.write(f"ROC AUC: {roc_auc:.4f}\n")
        f.write(f"ROC AUC (Dummy): {roc_auc_dummy:.4f}\n")
    