import networkx as nx
from degree_balaning_functions import back_and_forth_max_flow
np.random.seed(1234)


def main():
    pos_df = pd.read_csv(snakemake.input.set_pos, sep="\t", header=None)
    G_pos = nx.from_pandas_edgelist(pos_df, 0, 1, create_using=nx.DiGraph)

    G_comp = nx.DiGraph()
    G_comp.add_nodes_from(G_pos.nodes())

    for u in G_pos.nodes():
        for v in G_pos.nodes():
            if u != v and not G_pos.has_edge(u, v):
                G_comp.add_edge(u, v)
    
    selected_G, _ = back_and_forth_max_flow(G_pos, G_comp)
    nx.write_edgelist(
        selected_G[1], snakemake.output.train_neg, delimiter="\t", data=False
    )

main()