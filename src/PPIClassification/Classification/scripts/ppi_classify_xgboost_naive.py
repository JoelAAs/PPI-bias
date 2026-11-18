import sys
from classification_support import hyperparameter_tuned_model
import pandas as pd
from skopt.space import Integer
from xgboost import XGBClassifier
import joblib

# Same as ppi_classify_xgboost.py, but reads a merged bait/prey/interaction training set
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

RANDOM_STATE = snakemake.params.random_state
threads = snakemake.threads
directed = {"directional": True, "undirectional": False}[snakemake.wildcards.network_type]

df_train = pd.read_csv(snakemake.input.train, sep="\t", dtype={"bait": str, "prey": str})
df_pos = df_train.loc[df_train["interaction"], ["bait", "prey"]]
df_negative = df_train.loc[~df_train["interaction"], ["bait", "prey"]]

param_dist = [
    Integer(3, 10, name="max_depth"),
    Integer(100, 1000, name="n_estimators"),
    ]
static_parameters = dict(
    subsample=0.8
)

param_file = open(snakemake.output.params, "w")
model, parameters = hyperparameter_tuned_model(
    df_pos, df_negative,
    threads, directed,
    param_file,
    XGBClassifier,
    param_dist,
    static_parameters,
    snakemake.input.protein_embeddings,
    RANDOM_STATE,
    skip_hyper_parameter_tuning = True,
    default_parameters=False,
    n_iters=10)

joblib.dump(model, snakemake.output.saved_model)
param_file.close()
