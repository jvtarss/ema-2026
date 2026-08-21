import re
import pandas as pd

MASTER_EDGES_CSV = "master_edges.csv"
PPI_EDGES_CSV = "ppi_edges_harmonized_dedup.csv"
EDGES_OUTPUT_CSV = "network_edges_final.csv"
NODES_OUTPUT_CSV = "node_attributes.csv"

NOVEL_FAMILY_PATTERN = re.compile(r'^Egr-miRN\d')
KNOWN_PATTERN = re.compile(r'^Egr-miR\d')
NOVEL_LOCATION_PATTERN = re.compile(r'^(Chr\d+_\d+|scaffold_\d+_\d+)$')


def classify_mirna(node_id):
    if NOVEL_FAMILY_PATTERN.match(node_id):
        return "novel"
    if KNOWN_PATTERN.match(node_id):
        return "known"
    if NOVEL_LOCATION_PATTERN.match(node_id):
        return "novel"
    return None


def load_mirna_target_edges(path):
    df = pd.read_csv(path)
    mirna_edges = df[df["interaction_type"] == "miRNA-target"].copy()
    return mirna_edges[["source", "target", "score", "interaction_type", "inhibition_type"]]


def load_ppi_edges(path):
    df = pd.read_csv(path)
    return df[["source", "target", "score", "interaction_type", "inhibition_type"]]


def build_node_attributes(all_edges):
    all_nodes = sorted(set(all_edges["source"]) | set(all_edges["target"]))
    rows = []
    for node in all_nodes:
        status = classify_mirna(node)
        if status:
            rows.append({"node_id": node, "node_type": "miRNA", "mirna_status": status})
        else:
            rows.append({"node_id": node, "node_type": "protein", "mirna_status": "NA"})
    return pd.DataFrame(rows)


def main():
    mirna_edges = load_mirna_target_edges(MASTER_EDGES_CSV)
    ppi_edges = load_ppi_edges(PPI_EDGES_CSV)

    final_edges = pd.concat([mirna_edges, ppi_edges], ignore_index=True)
    final_edges = final_edges.drop_duplicates(subset=["source", "target", "interaction_type"])
    final_edges.to_csv(EDGES_OUTPUT_CSV, index=False)

    node_table = build_node_attributes(final_edges)
    node_table.to_csv(NODES_OUTPUT_CSV, index=False)

    print(len(final_edges), len(node_table))
    print(node_table["mirna_status"].value_counts())


if __name__ == "__main__":
    main()
