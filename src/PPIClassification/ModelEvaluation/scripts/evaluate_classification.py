import sys
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from evaluate_model import get_embedding_dict, centred_unit_embeddings, adjusted_roc


def read_pairs(file, predef):
    if predef:
        df = pd.read_csv(file, dtype={"protein1": str, "protein2": str})
        df = df.rename(columns={"protein1": "bait", "protein2": "prey"})
        df["interaction"] = df.pop("label").astype(bool)
    else:
        df = pd.read_csv(file, sep="\t", dtype={"bait": str, "prey": str})
    return df

def get_dataset(df, embedding_dict, flip_and_double):
    """
    Features from a merged bait/prey/interaction set, undirected pairs are scored in both orientations
    """
    pairs = df[["bait", "prey"]]
    y = df["interaction"].to_numpy(dtype=np.int8)
    if flip_and_double:
        flipped = pairs.rename(columns={"bait": "prey", "prey": "bait"})
        pairs = pd.concat([pairs, flipped], ignore_index=True)
        y = np.concatenate([y, y])
    X = np.hstack([
        np.vstack([embedding_dict[b] for b in pairs["bait"]]),
        np.vstack([embedding_dict[p] for p in pairs["prey"]])
    ]).astype(np.float32)
    return X, y, pairs


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    wc = snakemake.wildcards
    flip_and_double = {"directional": False, "undirectional": True}[wc.network_type]

    #states
    get_set_type = {"cand":"exp", "nonobs":"no"}
    predef = snakemake.params.predef
    if predef:
        train_no = get_set_type[snakemake.wildcards.state]
        test_no = get_set_type[snakemake.wildcards.test_state]
        seed = 0
    else:
        train_no = "no" if wc.no else "exp"
        test_no = "no" if wc.test_no else "exp"
        seed = wc.seed

    embedding_dict, _ = get_embedding_dict(snakemake.input.protein_embeddings)
    unit_emb = centred_unit_embeddings(embedding_dict)
    model = joblib.load(snakemake.input.saved_model)

    train_df = read_pairs(snakemake.input.train, predef)
    test_df = read_pairs(snakemake.input.test, predef)
    X_test, y_test, pairs_test = get_dataset(test_df, embedding_dict, flip_and_double)

    y_pred = model.predict_proba(X_test)[:, 1]
    y_pred_label = (y_pred >= 0.5).astype(np.int8)
    _, _, roc_auc_adj = adjusted_roc(y_test, y_pred, pairs_test, unit_emb)

    metrics = pd.DataFrame([{
        "network_type": wc.network_type,
        "positive": wc.pos_set,
        "negative": f"{wc.neg_set}-{train_no}",
        "seed": seed,
        "split": "predef" if predef else "naive",
        "esm_model": wc.esm_model,
        "classifier": wc.classifier,
        "roc_auc": roc_auc_score(y_test, y_pred),
        "roc_auc_adj": roc_auc_adj,
        "pr_auc": average_precision_score(y_test, y_pred),
        "acc_I": (y_pred_label[y_test == 1] == 1).mean(),
        "acc_NI": (y_pred_label[y_test == 0] == 0).mean(),
        "n_train_pos": int(train_df["interaction"].sum()),
        "n_train_neg": int((~train_df["interaction"]).sum()),
        "n_test_pos": int(test_df["interaction"].sum()),
        "n_test_neg": int((~test_df["interaction"]).sum()),
        "test_set": test_no
    }])
    metrics.to_csv(snakemake.output.metrics, sep="\t", index=False)
    pairs_test.assign(label=y_test, score=y_pred).to_csv(snakemake.output.pred, sep="\t", index=False)


if __name__ == "__main__":
    main()
