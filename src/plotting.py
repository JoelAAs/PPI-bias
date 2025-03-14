import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def max_observed_count(input_file, output_file):
    obs_max = dict()
    header = True
    with open(input_file, "r") as f:
        for line in f:
            if header:
                header = False
            else:
                _, _, max_interactions, observed_interactions = line.strip().split("\t")
                key = f"{max_interactions}:{observed_interactions}"
                if key in obs_max:
                    obs_max[key] += 1
                else:
                    obs_max[key] = 1

    with open(output_file, "w") as w:
        w.write(f"max_interactions\tobserved_interactions\tcount\n")
        for key in obs_max:
            count = obs_max[key]
            max_interactions, observed_interactions = [int(float(value)) for value in key.split(":")]
            w.write(f"{max_interactions}\t{observed_interactions}\t{count}\n")


observed_max_interaction_df = pd.read_csv("work_folder/observed_vs_possible.csv", sep = "\t")
max_interactions = np.max(observed_max_interaction_df["max_interactions"])
g = sns.JointGrid(
    x="max_interactions",
    y="observed_interactions",
    data=observed_max_interaction_df,
    xlim = (1, max_interactions),
    ylim = (1, max_interactions)
)

g.plot(sns.histplot, sns.histplot)
plt.show()