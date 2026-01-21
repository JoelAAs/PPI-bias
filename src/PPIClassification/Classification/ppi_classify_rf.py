import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import randint
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import ParameterSampler

global RANDOM_STATE


def get_dataset(pos_data_file, neg_data_file, embedding_dict, embed_length):
    df_pos = pd.read_csv(pos_data_file, sep="\t")
    df_negative = pd.read_csv(neg_data_file, sep="\t")
    n_samples = df_pos.shape[0] + df_negative.shape[0]
    X = np.zeros((n_samples, embed_length*2,), dtype=float)
    i = 0
    for c_df in [df_pos, df_negative]:
        for _, (bait_gene_name, prey_gene_name) in c_df.iterrows():
            X[i,:embed_length] = embedding_dict[bait_gene_name]
            X[i, embed_length:embed_length*2] = embedding_dict[prey_gene_name]
            i += 1
    y = np.zeros(n_samples, dtype=int)
    y[:df_pos.shape[0]] = 1
    return X, y




def get_embedding_dict(protein_embeddings_file):
    df_embed = pd.read_csv(protein_embeddings_file, sep="\t")
    embed_length = df_embed.shape[1] -1
    embedding_dict = {
        row.iloc[-1]: row.iloc[:-1].values  for i, row in df_embed.iterrows()
    }
    return embedding_dict, embed_length


def hyperparameter_tuned_model(X_train, y_train, X_test, y_test, threads):
    param_dist = {
        "n_estimators": randint(200, 100000),
        "max_depth": [None] + list(range(5, 50, 5)),
        "min_samples_split": randint(2, 20),
        "min_samples_leaf": randint(1, 10),
        "max_features": ["sqrt", "log2"]
    }

    best_score = -np.inf
    best_model = None
    best_params = None

    for params in ParameterSampler(param_dist, n_iter=30, random_state=RANDOM_STATE):
        model = RandomForestClassifier(
            **params,
            random_state=RANDOM_STATE,
            n_jobs=threads
        )

        model.fit(X_train, y_train)

        y_test_pred = model.predict(X_test)
        score = accuracy_score(y_test, y_test_pred)

        if score > best_score:
            best_score = score
            best_model = model
            best_params = params

    print("Best validation score:", best_score)
    print("Best params:", best_params)

    return best_model, score, params


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description="")
    # parser.add_argument("--train_ppi_data_pos", required=True, help="")
    # parser.add_argument("--train_ppi_data_neg", required=True, help="")
    # parser.add_argument("--test_ppi_data_pos", required=True, help="Path to output csv file")
    # parser.add_argument("--test_ppi_data_neg", required=True, help="Path to output csv file")
    # parser.add_argument("--validation_ppi_data_pos", required=True, help="Path to output csv file")
    # parser.add_argument("--validation_ppi_data_neg", required=True, help="Path to output csv file")
    #
    # parser.add_argument("--protein_embeddings", required=True, help="Path to output csv file")
    #
    # parser.add_argument("--threads", type=int, default=40, help="")
    # parser.add_argument("--randomstate", type=int, default=1234, help="")
    # args = parser.parse_args()
    #
    # RANDOM_STATE=args.randomstate
    # threads = args.threads
    # train_ppi_data_pos
    # train_ppi_data_neg
    # test_ppi_data_pos
    # test_ppi_data_neg
    # protein_embedings

    # embed_dict, n_embedding =  get_embedding_dict(args.protein_embeddings)
    #
    # X_train, y_train = get_dataset(
    #     args.train_ppi_data_pos,
    #     args.train_ppi_data_neg,
    #     embed_dict,
    #     n_embedding
    # )
    #
    # X_test, y_test = get_dataset(
    #     args.test_ppi_data_pos,
    #     args.test_ppi_data_neg,
    #     embed_dict,
    #     n_embedding
    # )
    #

    pn = "/per_gene"
    dataset="ms"
    threads = 40
    RANDOM_STATE=1234
    embed_dict, n_embedding =  get_embedding_dict(f"work_folder{pn}/embeddings/canonical_embedding.csv.gz")

    X_train, y_train = get_dataset(
        f"work_folder{pn}/subsets/train/balanced/{dataset}_pos.csv",
        f"work_folder{pn}/subsets/train/balanced/{dataset}_neg.csv",
        embed_dict,
        n_embedding
    )

    X_test, y_test = get_dataset(
        f"work_folder{pn}/subsets/test/balanced/{dataset}_pos.csv",
        f"work_folder{pn}/subsets/test/balanced/{dataset}_neg.csv",
        embed_dict,
        n_embedding
    )

    X_validate, y_validate = get_dataset(
        f"work_folder{pn}/subsets/test/balanced/{dataset}_pos.csv",
        f"work_folder{pn}/subsets/test/balanced/{dataset}_neg.csv",
        embed_dict,
        n_embedding
    )



    _, score, parameters = hyperparameter_tuned_model(X_train, y_train,X_test, y_test, threads)

    rfc = RandomForestClassifier(
        **parameters,
        n_jobs=threads)

    rfc.fit(
        np.vstack((X_train, X_test)),
        np.concatenate((y_train, y_test))
    )


    y_validate_pred = rfc.predict(X_validate)

    print("Final test accuracy:", accuracy_score(y_validate, y_validate_pred))
    print(classification_report(y_test, y_validate_pred))


