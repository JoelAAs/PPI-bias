import re

import joblib
import numpy as np
import pandas as pd

MODEL_PATH_RE = re.compile(r"_limit_(?P<neg_limit>[0-9.]+)_poslim_(?P<pos_limit>all|[0-9.]+)(?P<random>-random)?_model_")


def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] - 1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


def build_features(df, embedding_dict, embed_length, flip_and_double):
    baits = df.iloc[:, 0].to_numpy()
    prey = df.iloc[:, 1].to_numpy()
    n_samples = df.shape[0]

    if not flip_and_double:
        X = np.zeros((n_samples, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in baits]
        X[:, embed_length:] = [embedding_dict[p] for p in prey]
        pairs = df[["bait", "prey"]].reset_index(drop=True)
    else:
        X = np.zeros((n_samples * 2, embed_length * 2), dtype=np.float32)
        X[:, :embed_length] = [embedding_dict[b] for b in np.concatenate([baits, prey])]
        X[:, embed_length:] = [embedding_dict[p] for p in np.concatenate([prey, baits])]
        flipped = df.rename(columns={"bait": "prey", "prey": "bait"})[["bait", "prey"]]
        pairs = pd.concat([df[["bait", "prey"]], flipped], ignore_index=True)

    return X, pairs

if __name__ == "__main__":
    network_type = snakemake.wildcards.network_type
    if network_type == "undirectional":
        flip_and_double = True
    elif network_type == "directional":
        flip_and_double = False
    else:
        raise ValueError(f"{network_type} is an invalid network value")

    embedding_dict, embed_length = get_embedding_dict(snakemake.input.protein_embeddings)

    def parse_model_config(model_file):
        match = MODEL_PATH_RE.search(model_file)
        if not match:
            raise ValueError(f"Could not parse pos_limit/neg_limit/random from {model_file}")
        return match["pos_limit"], match["neg_limit"], match["random"] or ""


    models = {
        parse_model_config(model_file): joblib.load(model_file)
        for model_file in snakemake.input.models
    }

    datasets = [
        (snakemake.input.test_pos, "Interaction"),
        (snakemake.input.selected_negative, "Negative_HRNI"),
        (snakemake.input.no_selected_negative, "Negative_NO"),
    ]

    results = []
    for data_file, label in datasets:
        df = pd.read_csv(data_file, sep="\t", dtype={"bait": "string", "prey": "string"})
        X, pairs = build_features(df, embedding_dict, embed_length, flip_and_double)
        pairs["dataset"] = label
        for (pos_limit, neg_limit, random), model in models.items():
            train_negative = "no" if random == "-random" else "hrni"
            column = f"train_{train_negative}_poslim{pos_limit}_neglim{neg_limit}"
            pairs[column] = model.predict_proba(X)[:, 1]
        results.append(pairs)

    result_df = pd.concat(results, ignore_index=True)
    result_df.to_csv(snakemake.output.predictions, sep="\t", index=False)
