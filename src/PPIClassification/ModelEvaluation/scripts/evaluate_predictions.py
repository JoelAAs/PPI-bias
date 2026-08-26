import joblib
import numpy as np
import pandas as pd


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


def predict_dataset(data_file, dataset_label, embedding_dict, embed_length, flip_and_double, model_hrni, model_no):
    df = pd.read_csv(data_file, sep="\t", dtype={"bait": "string", "prey": "string"})
    X, pairs = build_features(df, embedding_dict, embed_length, flip_and_double)

    pairs["HRNI"] = model_hrni.predict_proba(X)[:, 1]
    pairs["NO"] = model_no.predict_proba(X)[:, 1]
    pairs["dataset"] = dataset_label
    return pairs


network_type = snakemake.wildcards.network_type
if network_type == "undirectional":
    flip_and_double = True
elif network_type == "directional":
    flip_and_double = False
else:
    raise ValueError(f"{network_type} is an invalid network value")

embedding_dict, embed_length = get_embedding_dict(snakemake.input.protein_embeddings)

model_hrni = joblib.load(snakemake.input.hrni_saved_model)
model_no = joblib.load(snakemake.input.no_saved_model)

datasets = [
    (snakemake.input.test_pos, "Interaction"),
    (snakemake.input.selected_negative, "Negative_HRNI"),
    (snakemake.input.non_obs, "Negative_NO"),
]

results = [
    predict_dataset(data_file, label, embedding_dict, embed_length, flip_and_double, model_hrni, model_no)
    for data_file, label in datasets
]

result_df = pd.concat(results, ignore_index=True)
result_df.to_csv(snakemake.output.predictions, sep="\t", index=False)
