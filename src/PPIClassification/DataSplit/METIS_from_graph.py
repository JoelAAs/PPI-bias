import argparse
import networkx as nx

def write_metis(graph_file, output_metis, gene_id_file, node_weights=False):
    G = nx.load_graphml(graph_file)

    sorted_nodes = list(enumerate(sorted(G.nodes())))
    int_mapping = {node: i + 1 for i, node in sorted_nodes}  # +1 for 1 index

    with open(gene_id_file, "w") as w:
        w.write("gene_name\tint_id\n")
        for gene, int_id in int_mapping.items():
            w.write(f"{gene}\t{int_id}\n")

    with open(output_metis, "w") as w:
        w.write(f'{G.number_of_nodes()} {G.number_of_edges()} 11\n')  # 11 for node and edge weights
        for i, node in sorted_nodes:
            if node_weights:
                line = f"{node["node_weight"]}  "
            else:
                line = ""
            line += " ".join([
                f'{int_mapping[edge[1]]} {edge[2]["edge_weight"]}' for
                edge in G.edges(node, data=True)])
            w.write(line + "\n")

if __name__ == '__main__':
    args = argparse.ArgumentParser(description="METIS from edge list")
    args.add_argument("--graph", "-g", required=True, help="edge list")
    args.add_argument("--output_metis", "-g", required=True, help="mentis output file")
    args.add_argument("--output_int_id", "-o", required=True, help="output file")
    args.add_argument("--node_weights","-w", default=False, action="store_true")

    graph = args.graph
    output_metis = args.output_metis
    output_int_id = args.output_int_id
    node_weights = args.node_weights
    write_metis(graph, output_metis, output_int_id, node_weights)
