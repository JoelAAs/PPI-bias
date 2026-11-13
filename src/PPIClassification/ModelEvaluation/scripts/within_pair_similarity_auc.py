import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def load_embeddings(path):
    df = pd.read_csv(path, sep="\t").set_index("id_col")
    E = df.to_numpy(np.float64)
    E = E - E.mean(axis=0)  # centre: remove the direction shared by all proteins
    E /= np.linalg.norm(E, axis=1, keepdims=True)  # unit length, so cosine = dot product
    return dict(zip(df.index, E))


def auc(pos_sim, neg_sim):
    y = np.r_[np.ones(len(pos_sim)), np.zeros(len(neg_sim))]
    return roc_auc_score(y, np.r_[pos_sim, neg_sim])


def cosim_emb(path, embedding_dict):
    df = pd.read_csv(path, sep="\t", usecols=["bait", "prey"])
    return np.array([np.dot(embedding_dict[b], embedding_dict[p]) for b, p in zip(df["bait"], df["prey"])])


def main():
    embedding_dict = load_embeddings(snakemake.input.protein_embeddings)

    train_cossim_pos = cosim_emb(snakemake.input.train_pos, embedding_dict)
    train_cossim_hrni = cosim_emb(snakemake.input.train_hrni, embedding_dict)
    train_cossim_no = cosim_emb(snakemake.input.train_no, embedding_dict)

    test_cossim_pos = cosim_emb(snakemake.input.test_pos, embedding_dict)
    test_cossim_hrni = cosim_emb(snakemake.input.test_hrni, embedding_dict)
    test_cossim_no = cosim_emb(snakemake.input.test_no, embedding_dict)

    row = {
        "source": snakemake.wildcards.source,
        "network_type": snakemake.wildcards.network_type,
        "esm_model": snakemake.wildcards.esm_model,
        "dataset": snakemake.wildcards.dataset,
        "pos_limit": snakemake.wildcards.pos_limit,
        "neg_limit": snakemake.wildcards.neg_limit,
        "train_hrni": auc(train_cossim_pos, train_cossim_hrni),
        "train_no": auc(train_cossim_pos, train_cossim_no),
        "test_hrni": auc(test_cossim_pos, test_cossim_hrni),
        "test_no": auc(test_cossim_pos, test_cossim_no),
    }
    pd.DataFrame([row]).to_csv(snakemake.output.auc, sep="\t", index=False)


main()
