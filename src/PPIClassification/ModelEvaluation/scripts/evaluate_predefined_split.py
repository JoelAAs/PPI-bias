import sys
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from evaluate_model import get_embedding_dict, centred_unit_embeddings, adjusted_roc
from evaluate_naive_split import get_dataset

# Predefined (ppi-split pipeline) state -> naive-split naming: cand = experimental negatives, nonobs = NO negatives
STATE_SUFFIX = {"cand": "", "nonobs": "-no"}


def read_split(path):
    df = pd.read_csv(path, dtype={"protein1": str, "protein2": str})
    df = df.rename(columns={"protein1": "bait", "protein2": "prey"})
    df["interaction"] = df.pop("label").astype(bool)
    return df


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    wc = snakemake.wildcards
    flip_and_double = {"directional": False, "undirectional": True}[wc.network_type]

    embedding_dict, _ = get_embedding_dict(snakemake.input.protein_embeddings)
    unit_emb = centred_unit_embeddings(embedding_dict)
    model = joblib.load(snakemake.input.saved_model)

    train_df = read_split(snakemake.input.train)
    test_df = read_split(snakemake.input.test)
    X_test, y_test, pairs_test = get_dataset(test_df, embedding_dict, flip_and_double)

    y_pred = model.predict_proba(X_test)[:, 1]
    y_pred_label = (y_pred >= 0.5).astype(np.int8)
    _, _, roc_auc_adj = adjusted_roc(y_test, y_pred, pairs_test, unit_emb)

    metrics = pd.DataFrame([{
        "network_type": wc.network_type,
        "positive": wc.pos_set,
        "negative": f"{wc.neg_set}{STATE_SUFFIX[wc.state]}",
        "seed": "predefined",
        "esm_model": wc.esm_model,
        "classifier": "xgboost",
        "roc_auc": roc_auc_score(y_test, y_pred),
        "roc_auc_adj": roc_auc_adj,
        "pr_auc": average_precision_score(y_test, y_pred),
        "acc_I": (y_pred_label[y_test == 1] == 1).mean(),
        "acc_NI": (y_pred_label[y_test == 0] == 0).mean(),
        "n_train_pos": int(train_df["interaction"].sum()),
        "n_train_neg": int((~train_df["interaction"]).sum()),
        "n_test_pos": int(test_df["interaction"].sum()),
        "n_test_neg": int((~test_df["interaction"]).sum()),
        "test_set": "no" if wc.test_state == "nonobs" else "exp",
        "train_state": wc.state
    }])
    metrics.to_csv(snakemake.output.metrics, sep="\t", index=False)
    pairs_test.assign(label=y_test, score=y_pred).to_csv(snakemake.output.pred, sep="\t", index=False)


if __name__ == "__main__":
    main()
