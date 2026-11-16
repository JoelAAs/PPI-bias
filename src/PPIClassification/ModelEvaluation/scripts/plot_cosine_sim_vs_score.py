import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def load_unit_embeddings(path):
    df = pd.read_csv(path, sep="\t").set_index("id_col")
    embeddings = df.to_numpy(np.float64)
    embeddings = embeddings - embeddings.mean(axis=0)  # centre: remove the direction shared by all proteins
    embeddings /= np.linalg.norm(embeddings, axis=1, keepdims=True)  # unit length, so cosine = dot product
    return df.index.to_numpy(), embeddings.astype(np.float32)


def cosim_emb(bait, prey, ids, embeddings):
    sim = embeddings @ embeddings.T  # all protein-protein cosine similarities (~20k x 20k)
    id_set = pa.array(ids)
    b = pc.index_in(bait, value_set=id_set).to_numpy()
    p = pc.index_in(prey, value_set=id_set).to_numpy()
    return sim[b, p]


def get_values(score, cos_sim):
    q25, q50, q75 = np.percentile(cos_sim, [25, 50, 75])
    return {
        "range_median": np.median(score),
        "cos_sim_median": q50,
        "n_pairs": len(score),
        "q25": q25,
        "q75": q75,
    }


def thin(table, score_col, n_points=5000):
    # too many rows to draw: keep rows evenly spaced in rank and evenly spaced in score (table sorted ascending)
    score = table[score_col].to_numpy()
    idx = np.unique(np.concatenate([
        np.linspace(0, len(score) - 1, n_points).astype(int),
        np.searchsorted(score, np.linspace(score[0], score[-1], n_points)),
    ]))
    return table.iloc[idx]


def plot_stats(curves, negative, score_col, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))

    for col, color, label in [
        ("cum_mean_max", "#2a6fdb", f"pairs with {score_col} ≤ x"),
        ("cum_mean_min", "#d9480f", f"pairs with {score_col} ≥ x"),
    ]:
        sd = np.sqrt(curves[f"{col}_var"])
        ax.fill_between(curves[score_col], curves[col] - sd, curves[col] + sd, color=color, alpha=0.15, linewidth=0)
        ax.plot(curves[score_col], curves[col], color=color, linewidth=1.2, label=label)

    ax.errorbar(
        negative["x"], negative["mean"], yerr=np.sqrt(negative["var"]),
        fmt="o", color="black", markersize=5, capsize=3, zorder=3,
        label=f"negatives (n = {negative['n']:,})",
    )
    ax.axhline(0, color="0.6", linewidth=0.8, linestyle="--", zorder=0)

    ax.set_xlabel(score_col)
    ax.set_ylabel("Cumulative mean cosine similarity (± SD)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    score_col = snakemake.params.score_col
    ids, embeddings = load_unit_embeddings(snakemake.input.protein_embeddings)

    bait_col, prey_col = snakemake.params.bait_col, snakemake.params.prey_col
    id_columns = [bait_col, prey_col, score_col]
    if snakemake.wildcards.source == "pod":
        id_columns += ["n_tested", "n_observed"]
    else:
        id_columns.append("posterior")
        
    table = pq.read_table(snakemake.input.dataset, columns=id_columns)
    cos_sim = cosim_emb(table[bait_col], table[prey_col], ids, embeddings)
    table = table.select(id_columns[2:]).to_pandas()  # drop the id strings before going to pandas
    table["cos_emb"] = cos_sim
    
    neg_limit = snakemake.params.neg_limit
    if snakemake.wildcards.source == "pod":
        is_neg = (table["n_observed"] == 0) & (table["n_tested"] >= neg_limit)
    elif snakemake.wildcards.source == "fdr":
        is_neg = table["posterior"] <= neg_limit
    neg_table, table = table[is_neg], table[~is_neg]  # negatives -> one dot, the rest -> cumulative curves
    
    negative = {
        "x": neg_table[score_col].median(),
        "mean": neg_table["cos_emb"].mean(),
        "var": neg_table["cos_emb"].var(),
        "n": len(neg_table),
    }
    
    table = table.sort_values(score_col)
    table["cum_mean_max"] = table["cos_emb"].expanding().mean()
    table["cum_mean_max_var"] = table["cos_emb"].expanding().var()
    table = table.sort_values(score_col, ascending=False)
    table["cum_mean_min"] = table["cos_emb"].expanding().mean()
    table["cum_mean_min_var"] = table["cos_emb"].expanding().var()

    curves = thin(table.iloc[::-1], score_col)  # back to ascending score
    plot_stats(curves, negative, score_col, snakemake.output.emb_sim_vs_score)
    

main()
