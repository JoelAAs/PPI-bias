import pandas as pd
import networkx as nx

def get_spearman(G_pos, G_neg):
    bait_degree = [dict(G_pos.out_degree()), dict(G_neg.out_degree())]
    prey_degree = [dict(G_pos.in_degree()), dict(G_neg.in_degree())]
    
    baits = set(bait_degree[0].keys()) | set(bait_degree[1].keys())
    df_baits = pd.DataFrame({"node": list(baits)})
    df_baits["type"] = "bait"
    df_baits["pos_degree"] = df_baits["node"].apply(lambda x: bait_degree[0].get(x,0))
    df_baits["neg_degree"] = df_baits["node"].apply(lambda x: bait_degree[1].get(x,0))
    
    
    prey =  set(prey_degree[0].keys()) | set(prey_degree[1].keys())
    df_prey = pd.DataFrame({"node": list(prey)})
    df_prey["type"] = "prey"
    df_prey["pos_degree"] = df_prey["node"].apply(lambda x: prey_degree[0].get(x,0))
    df_prey["neg_degree"] = df_prey["node"].apply(lambda x: prey_degree[1].get(x,0))
    
    df_degree = pd.concat([df_baits, df_prey], axis=0)

    spearman = df_degree[["pos_degree","neg_degree"]].corr(method="spearman").iloc[0,1]
    
    return spearman

def get_bait_prey_div(G_pos, G_neg):
    bait_degree = [dict(G_pos.out_degree()), dict(G_neg.out_degree())]
    prey_degree = [dict(G_pos.in_degree()), dict(G_neg.in_degree())]

    df_delta = pd.DataFrame({"node": list(set(G_pos.nodes()) | set(G_neg.nodes()))})
    df_delta["pos_bait"] = df_delta["node"].apply(lambda x: bait_degree[0].get(x,0))
    df_delta["neg_bait"] = df_delta["node"].apply(lambda x: bait_degree[1].get(x,0))
    df_delta["pos_prey"] = df_delta["node"].apply(lambda x: prey_degree[0].get(x,0))
    df_delta["neg_prey"] = df_delta["node"].apply(lambda x: prey_degree[1].get(x,0))

    df_delta["delta_bait"] = (df_delta["pos_bait"] - df_delta["neg_bait"]).abs()
    df_delta["delta_prey"] = (df_delta["pos_prey"] - df_delta["neg_prey"]).abs()

    return df_delta["delta_bait"].sum(), df_delta["delta_prey"].sum()


def read_graphs(input_files):
    return [
        nx.from_pandas_edgelist(pd.read_csv(f, sep="\t",header=None), 0, 1, create_using=nx.DiGraph)
        for f in input_files
        ]

def main():
    source_graphs = read_graphs([snakemake.input.set_pos, snakemake.input.set_neg])
    train_graphs = read_graphs([snakemake.input.train_pos, snakemake.input.train_neg])
    validation_graphs = read_graphs([snakemake.input.validation_pos, snakemake.input.validation_neg])
    test_graphs = read_graphs([snakemake.input.test_pos, snakemake.input.test_neg])
    
    original_edges = source_graphs[0].number_of_edges()
    with open(snakemake.output.edge_statistics, "w") as w:
        w.write("dataset\tsplit\tspearman\tdivergence\tn_edges\tfraction_total\n")
        for graphs, split_name in zip([
            source_graphs,
            train_graphs,
            validation_graphs,
            test_graphs
        ], [
            "source",
            "train",
            "validation",
            "test"
        ]):
            dataset = f"{snakemake.wildcards.dataset}_{snakemake.wildcards.neg_limit}_{snakemake.wildcards.pos_limit}"
            spearman = get_spearman(*graphs)
            div_degrees = sum(get_bait_prey_div(*graphs))
            n_edges = graphs[0].number_of_edges() # assume n_pos ~ n_neg where it matters
            fraction_total = n_edges/original_edges
            w.write(f"{dataset}\t{split_name}\t{spearman}\t{div_degrees}\t{n_edges}\t{fraction_total}\n")

main()