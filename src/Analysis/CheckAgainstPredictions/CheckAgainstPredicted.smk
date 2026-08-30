import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
from collections import defaultdict



rule get_zhang_heteromer:
    """
    https://datadryad.org/dataset/doi:10.5061/dryad.15dv41p84
    needs final_predictions.tar.gz. download manually to work_folder/data/predicted_interactions and un-tar it
    """
    output:
        "work_folder/data/predicted_interactions/final_predictions/final_predictions_80.tsv",
        "work_folder/data/predicted_interactions/final_predictions/final_predictions_90.tsv",
        "work_folder/data/predicted_interactions/RF2-PPI_scores",
        "work_folder/data/predicted_interactions/AF_scores",
        "work_folder/data/predicted_interactions/DCA_scores"
    run:
        
        raise IOError("You need to download these files yourself, cant wget")


rule join_top_to_negatome:
    # prediction columns:
    # RFprob — RF2-PPI contact probability (the fast screening network).
    # AFprob — AlphaFold2 contact probability, model 3, run on their omicMSAs.
    # CFprob — ColabFold pipeline using the AF2 network.
    # AFMprob — ColabFold using the AlphaFold-Multimer network.
    input:
        interaction_prediction = "work_folder/data/predicted_interactions/final_predictions/final_predictions_80.tsv",
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq" # Obs undirectional
    output:
        joined = "work_folder/analysis/fpr_comparission/predicted/Zhang_{dataset}_joined.pq"
    run:
        df_af_predictions = pd.read_csv(input.interaction_prediction, sep="\t", comment="#")
        df_pod = pd.read_parquet(input.pod)

        df_pod["pair_id"] = df[["uniprot_id_bait","uniprot_id_prey"]].min(axis=1) + ":" + df[["a","b"]].max(axis=1)
        df_af_predictions["pair_id"] = df_af_predictions.apply(lambda row: ":".join(sorted([row["Protein1"], row["Protein2"]])), axis=1)

        df_merged = df_pod.merge(df_af_predictions, on="pair_id", how="inner")
        df_merged.to_parquet(output.joined)


rule sort_score:
    input: 
        scores = "work_folder/data/predicted_interactions/{model}_scores", # AF, RF2, DCA
    output:
        protein_sorted = temp("work_folder/data/predicted_interactions/{model}_scores_protein_sorted") # AF, RF2, DCA
    run:
        with open(output.protein_sorted, "w") as w, open(input.scores, "r") as f:
            for line in f:
                if line[0] == "#" or line.startswith("Pair"):
                    continue
                else:
                    cols = line.split("\t")
                    prot_id = "_".join(sorted(cols[0].split("_")))
                    w.write(f"{prot_id}\t{'\t'.join(cols[1:])}")

rule line_sort_score:
    input:
        protein_sorted = "work_folder/data/predicted_interactions/{model}_scores_protein_sorted"
    output:
        line_sorted = temp("work_folder/data/predicted_interactions/{model}_scores_sorted")
    threads:
        10
    shell:
        """
        if [ "{wildcards.model}" != "DCA" ]; then
            echo -e 'Pair\tinteraction_probability\tSource' > {output.line_sorted}
        else
            echo -e 'Pair\tinteraction_probability' > {output.line_sorted}
        fi
        LC_ALL=C sort --parallel={threads} -t$'\t' -k1,1 {input.protein_sorted} >> {output.line_sorted}
        """


rule join_predictions:
    input:
        scores = expand("work_folder/data/predicted_interactions/{model}_scores_sorted",
                        model=prediction_models),
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq"
    output:
        joined_scores = "work_folder/data/predicted_interactions/joined/{dataset}_scores.pq"
    resources:
        mem_gb = 50
    log:
        "logs/data/predicted_interactions/joined/{dataset}_scores.log"
    run:
        with open(log[0], "w") as fh:
            def note(m): fh.write(m + "\n"); fh.flush()

            df = pd.read_parquet(input.pod, columns=["uniprot_id_bait", "uniprot_id_prey",
                                                     "n_tested", "n_observed"])
            
            ids = df[["uniprot_id_bait", "uniprot_id_prey"]]
            df["pair_id"] = ids.min(axis=1) + "_" + ids.max(axis=1)
            df = df.sort_values("pair_id", ignore_index=True)
            pairs = df["pair_id"].tolist()
            n = len(pairs)
            note(f"Line in POD file: {n:,}")

            for path, m in zip(input.scores, prediction_models):
                scores = np.full(n, np.nan, dtype=np.float32)
                subset = np.empty(n, dtype=object) if m != "DCA" else None
                matched = unmatched = 0
                i = 0

                with open(path) as f:
                    for line in f:
                        cols = line.rstrip("\n").split("\t")
                        pid = cols[0]
                        if pid == "Pair":
                            continue
                        while i < n and pairs[i] < pid:
                            i += 1
                        if i >= n:
                            unmatched += 1
                            continue
                        if pairs[i] == pid:
                            scores[i] = float(cols[1])
                            if subset is not None:
                                subset[i] = cols[2]
                            matched += 1
                            i += 1
                        else:
                            unmatched += 1

                df[f"{m}_score"] = scores
                if subset is not None:
                    df[f"{m}_subset"] = subset
                note(f"{m}: matched {matched:,}, unmatched {unmatched:,}")

            df.to_parquet(output.joined_scores)
            note(f"joined: {len(df):,}")

rule survival_curves:
    input:
        joined = "work_folder/data/predicted_interactions/joined/{dataset}_scores.pq"
    output:
        curves = "work_folder/analysis/predicted_interactions/survival/{dataset}_survival.tsv",
        counts = "work_folder/analysis/predicted_interactions/survival/{dataset}_counts.tsv"
    params:
        n_thresholds = [1, 3, 5]
    run:
        DCA_EPS = 1e-6

        df = pd.read_parquet(input.joined)
        df = df[df["DCA_score"].notna()].copy()
        df["DCA_score"] = np.log(df["DCA_score"].clip(lower=DCA_EPS))

        is_neg = pd.Series(False, index=df.index)
        for m in ("RF2-PPI", "AF"):
            col = f"{m}_subset"
            if col in df.columns:
                is_neg |= df[col].eq("NEG")

        tested   = df["n_tested"].fillna(0)
        observed = df["n_observed"].fillna(-1)

        sets = {"random_NEG": df[is_neg]}
        for k in params.n_thresholds:
            sets[f"HRNI_n{k}"] = df[(observed == 0) & (tested >= k)]
        sets["HRI"] = df[observed > 0]

        grids = {
            "DCA": np.quantile(df["DCA_score"], np.linspace(0, 1, 201)),
            "RF2-PPI": np.linspace(0, 1, 201),
            "AF":  np.linspace(0, 1, 201),
        }

        rows, counts = [], []
        for set_name, sub in sets.items():
            n_elig = len(sub)
            for model, grid in grids.items():
                col = f"{model}_score"
                n_scored = sub[col].notna().sum()
                counts.append(dict(set=set_name, model=model,
                                   n_eligible=n_elig, n_scored=n_scored,
                                   coverage=n_scored / n_elig if n_elig else np.nan))
                uncond = sub[col].fillna(-np.inf).to_numpy()
                cond   = sub[col].dropna().to_numpy()
                for t in grid:
                    rows.append(dict(
                        set=set_name, model=model, threshold=float(t),
                        frac_above_uncond=float((uncond > t).mean()) if n_elig else np.nan,
                        frac_above_cond=float((cond > t).mean()) if n_scored else np.nan,
                        n_above=int((uncond > t).sum()),
                    ))

        pd.DataFrame(rows).to_csv(output.curves, sep="\t", index=False)
        pd.DataFrame(counts).to_csv(output.counts, sep="\t", index=False)


rule plot_survival:
    input:
        curves = "work_folder/analysis/predicted_interactions/survival/{dataset}_survival.tsv",
        counts = "work_folder/analysis/predicted_interactions/survival/{dataset}_counts.tsv"
    output:
        png = "work_folder/analysis/predicted_interactions/survival/plot/{dataset}_survival.png"
    params:
        # score thresholds to mark. RF2 0.3 is the stage-2 gate; the de novo
        # calling cutoffs sit near 0.99 (Fig. 3C/3D). DCA has no published
        # cutoff, so nothing is marked there.
        cutoffs = {"RF2-PPI": [0.3, 0.99], "AF": [0.99]},
        denom = "uncond"   # "uncond" for the paper, "cond" for the funnel-conditional version
    run:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        curves = pd.read_csv(input.curves, sep="\t")
        counts = pd.read_csv(input.counts, sep="\t")

        ycol = f"frac_above_{params.denom}"
        models = ["DCA", "RF2-PPI", "AF"]

        style = {
            "HRI":        dict(color="#2ca02c", ls="-",  lw=1.8, label="Interaction (any)"),
            "random_NEG": dict(color="#1f77b4", ls="--", lw=1.8, label=" Pre-selected Random negatives"),
            "HRNI_n1":    dict(color="#ff7f0e", ls="-",  lw=1.2, alpha=0.55, label="HRNI, $N\\geq1$"),
            "HRNI_n3":    dict(color="#ff7f0e", ls="-",  lw=1.5, alpha=0.78, label="HRNI, $N\\geq3$"),
            "HRNI_n5":    dict(color="#ff7f0e", ls="-",  lw=1.9, alpha=1.0,  label="HRNI, $N\\geq5$"),
        }

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

        for ax, model in zip(axes, models):
            sub_m = curves[curves["model"] == model]

            for set_name, st in style.items():
                s = sub_m[sub_m["set"] == set_name].sort_values("threshold")
                if s.empty:
                    continue
                y = s[ycol].to_numpy(dtype=float)
                # 0 cannot be drawn on a log axis; clip to the resolution
                # floor of that set so "zero" reads as "below 1/n", not as absent
                n = counts.loc[
                    (counts["set"] == set_name) & (counts["model"] == model),
                    "n_eligible"
                ]
                floor = 1.0 / n.iloc[0] if len(n) and n.iloc[0] else 1e-6
                ax.plot(s["threshold"], np.clip(y, floor, None), **st)

                # mark the floor so the reader knows where measurement stops
                ax.axhline(floor, color=st["color"], lw=0.5, ls=":", alpha=0.35)

            for t in params.cutoffs.get(model, []):
                ax.axvline(t, color="grey", lw=0.8, ls="-.", alpha=0.7)
                ax.text(t, 1.02, f"{t:g}", transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=7, color="grey")

            ax.set_yscale("log")
            ax.set_ylim(top=1.5)
            ax.set_title(model)
            ax.set_xlabel("log(DCA) score" if model == "DCA" else "Interaction probability")
            ax.grid(alpha=0.25, which="both", lw=0.4)

        axes[0].set_ylabel(
            "Fraction of pairs above threshold"
            + ("" if params.denom == "uncond" else "\n(among pairs scored by this model)")
        )

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False,
                   bbox_to_anchor=(0.5, -0.06))

        fig.savefig(output.png, dpi=300, bbox_inches="tight")
        plt.close(fig)