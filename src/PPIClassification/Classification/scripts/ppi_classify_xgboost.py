import sys
from classification_support import hyperparameter_tuned_model
import pandas as pd
from skopt.space import Integer, Real
from xgboost import XGBClassifier
import joblib

# keep the previous "> {log} 2>&1" behaviour
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

RANDOM_STATE = snakemake.params.random_state
threads = snakemake.threads
network_type = snakemake.wildcards.network_type
if network_type == "undirectional":
    directed = False
elif network_type == "directional":
    directed = True
else:
    raise ValueError(f"{network_type} is an invalid network value")

df_pos = pd.read_csv(snakemake.input.train_pos, sep="\t")[["bait", "prey"]]
df_negative = pd.read_csv(snakemake.input.train_neg, sep="\t")[["bait", "prey"]]

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
