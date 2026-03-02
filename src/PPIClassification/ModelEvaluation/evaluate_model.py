import pandas as pd 
import numpy as np 
import argparse
import joblib
from sklearn.metrics import precision_recall_curve, auc, roc_curve, log_loss
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

def generate_and_plot_performace(model, X_test, y_test, pr_png, neg_pr_png, roc_png, ce_png):
    y_pred = model.predict_proba(X_test)[:, 1]
    precision, recall, _ = precision_recall_curve(y_test, y_pred)
    pr_auc = auc(recall, precision)
    random_chance = np.mean(y_test)

    base_dist_pr, auc_base_dist_pr = get_baseline_performance(y_pred, y_test)
    print("plotting PR curve", flush=True)
    get_base_line_plot(
        (recall, precision), pr_auc, base_dist_pr, auc_base_dist_pr, pr_png, random_chance, eval_method_name="PR"
    )
    base_dist_roc, auc_base_dist_roc = get_baseline_performance(y_pred, y_test, eval_method=roc_curve)
    fpr, tpr, _ = roc_curve(y_test, y_pred)
    roc_auc = auc(fpr, tpr)
    print("plotting ROC curve", flush=True)
    get_base_line_plot(
        (fpr, tpr), roc_auc, base_dist_roc, auc_base_dist_roc, roc_png, None, eval_method_name="ROC"
    )

    y_pred_neg = 1 - y_pred
    precision_neg, recall_neg, _ = precision_recall_curve(1-y_test, y_pred_neg)
    pr_auc_neg = auc(recall_neg, precision_neg)
    random_chance_neg = np.mean(1-y_test)
    base_dist_pr_neg, auc_base_dist_pr_neg = get_baseline_performance(y_pred_neg, 1-y_test)
    print("plotting neg PR curve", flush=True)
    get_base_line_plot(
        (recall_neg, precision_neg), pr_auc_neg, base_dist_pr_neg, auc_base_dist_pr_neg, neg_pr_png, random_chance_neg, eval_method_name="PR"
    )

    dist_ce, obs_ce = get_baseline_log_loss(y_pred, y_test)
    get_cross_entropy_plot(dist_ce, obs_ce, ce_png)

    return pr_auc, pr_auc_neg, roc_auc, obs_ce, np.mean(auc_base_dist_pr), np.mean(auc_base_dist_pr_neg), np.mean(auc_base_dist_roc), np.mean(dist_ce)


def get_baseline_log_loss(y_pred, y_test, n = 1000):
    baseline_logloss = np.zeros(n, dtype=float)
    for i in range(n):
        y_pred_permut = np.random.permutation(y_pred)        
        baseline_logloss[i] = log_loss(y_test, y_pred_permut)
    
    return baseline_logloss, log_loss(y_test, y_pred)


def get_cross_entropy_plot(baseline_ce, obs_ce, ce_output):
    n_permutations = len(baseline_ce)

    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    axes[1].hist(baseline_ce, bins=int(math.sqrt(n_permutations)), alpha=0.7, label='Baseline CE Distribution')
    axes[1].axvline(obs_ce, color='blue', linestyle='--', label=f'Observed CE: {obs_ce:.4f}')
    axes[1].set_xlabel('Cross entropy')
    axes[1].set_title('Cross entropy Distribution Baseline')
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(ce_output)
    plt.close()
     


def get_baseline_performance(y_pred, y_test, eval_method=precision_recall_curve, n=1000):
    n_thresholds = len(set(y_pred))+1
    base_dist = np.zeros((n_thresholds*n, 3))
    auc_base_dist = []
    k = 0
    for i in range(n):
        y_pred_dummy_permut = np.random.permutation(y_pred)
        xi, yi, _ = eval_method(y_test, y_pred_dummy_permut)
        if eval_method == precision_recall_curve:
                x, y =  yi, xi # swap for PR curve to have recall on x-axis and precision on y-axis
        else:
                x, y = xi, yi
        sorted_idx = np.argsort(x)
        len_x = len(x)
        base_dist[k:k+len_x, :3] =np.column_stack([x, y, [i]*len_x])
        pr_auc = auc(x[sorted_idx], y[sorted_idx])
        auc_base_dist.append(pr_auc)
        k += len_x
    base_dist = base_dist[:k, :]
    return base_dist, auc_base_dist


def get_base_line_plot(obs_performance, obs_auc, base_dist, auc_base_dist, output_png, random_chance, eval_method_name="PR"):
    n_permutations = len(auc_base_dist)
    # Either precision-recall or FDR-TPR
    random_n_index = np.random.choice(range(base_dist.shape[0]), size=min(1000, base_dist.shape[0]), replace=False)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    axes[0].scatter(base_dist[random_n_index, 0], base_dist[random_n_index, 1], alpha=0.3, label="Permuted Performance", color='orange')
    axes[0].plot(obs_performance[0], obs_performance[1], color='blue', label='Observed Performance')
    if eval_method_name == "PR":
        axes[0].plot([0, 1], [random_chance, random_chance], linestyle='--', label='Prevalence')
    axes[0].set_xlabel('Recall' if eval_method_name == "PR" else 'FPR')
    axes[0].set_ylabel('Precision' if eval_method_name == "PR" else 'TPR')
    axes[0].set_title(f'{eval_method_name} Curve with Baseline')
    axes[0].legend()

    # AUC distribution
    axes[1].hist(auc_base_dist, bins=int(math.sqrt(n_permutations)), alpha=0.7, label='Baseline AUC Distribution')
    axes[1].axvline(obs_auc, color='blue', linestyle='--', label=f'Observed AUC: {obs_auc:.4f}')
    axes[1].set_xlabel(f'{eval_method_name} AUC')
    axes[1].set_title(f'{eval_method_name} AUC Distribution with Baseline')
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()

def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length):
    if pos_data_file[-3:] == "csv":
        df_pos = pd.read_csv(pos_data_file, sep="\t", usecols=[0, 1], header=None)
        df_negative = pd.read_csv(neg_data_file, sep="\t", usecols=[0, 1], header=None, comment="#")
    else:
        df_pos = pd.read_parquet(pos_data_file).iloc[:,:2]
        df_negative = pd.read_parquet(neg_data_file).iloc[:,:2]
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate a trained model on a test set")
    parser.add_argument("--pos_data_file", type=str, required=True, help="Path to the positive samples file")
    parser.add_argument("--neg_data_file", type=str, required=True, help="Path to the negative samples file")
    parser.add_argument("--protein_embeddings_file", type=str, required=True, help="Path to the protein embeddings file")
    parser.add_argument("--model_file", type=str, required=True, help="Path to the saved model file")
    parser.add_argument("--output_file", type=str, help="File to save evaluation results")   
    parser.add_argument("--plot_pr_png", type=str)   
    parser.add_argument("--plot_neg_pr_png", type=str)   
    parser.add_argument("--plot_roc_png", type=str)   
    parser.add_argument("--plot_ce_png", type=str)   


    args = parser.parse_args()


    embedding_dict, embed_length = get_embedding_dict(args.protein_embeddings_file)
    X_test, y_test = get_dataset(args.pos_data_file, args.neg_data_file, embedding_dict, embed_length)
    model = joblib.load(args.model_file)

    pr_auc, pr_auc_neg, roc_auc, obs_ce, base_pr_auc, base_pr_auc_neg, base_roc_auc, base_ce = generate_and_plot_performace(
        model, X_test, y_test, args.plot_pr_png,  args.plot_neg_pr_png,  args.plot_roc_png, args.plot_ce_png)
    
    with open(args.output_file, "w") as f:
        f.write(f"PR AUC: {pr_auc:.4f}\n")
        f.write(f"PR AUC (baseline): {base_pr_auc:.4f}\n")
        f.write(f"PR NEG AUC: {pr_auc_neg:.4f}\n")
        f.write(f"PR NEG AUC (baseline): {base_pr_auc_neg:.4f}\n")
        f.write(f"ROC AUC: {roc_auc:.4f}\n")
        f.write(f"ROC AUC (baseline): {base_roc_auc:.4f}\n")
        f.write(f"CE: {obs_ce:.4f}\n")
        f.write(f"CE (baseline): {base_ce:.4f}\n")
        f.write(f"Samples (pos/neg): {sum(y_test)} / {len(y_test)-sum(y_test)}")