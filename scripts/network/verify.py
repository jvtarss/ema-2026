import pandas as pd
import networkx as nx

EDGES_CSV = "network_edges_final.csv"
TOP_N_DEGREE = 10


def main():
    df = pd.read_csv(EDGES_CSV)
    G = nx.from_pandas_edgelist(df, source="source", target="target", edge_attr=True, create_using=nx.Graph())

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_mirna_target = (df["interaction_type"] == "miRNA-target").sum()
    n_ppi = (df["interaction_type"] == "ppi-string").sum()
    avg_degree = sum(dict(G.degree()).values()) / n_nodes if n_nodes else 0
    avg_clustering = nx.average_clustering(G)
    degree_sorted = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)

    print(n_nodes, n_edges, n_mirna_target, n_ppi, round(avg_degree, 3), round(avg_clustering, 3))
    for node, degree in degree_sorted[:TOP_N_DEGREE]:
        print(node, degree)


if __name__ == "__main__":
    main()
