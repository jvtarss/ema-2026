import pandas as pd

NODES_CSV = "node_attributes.csv"


def main():
    node_table = pd.read_csv(NODES_CSV)
    node_table["color_group"] = node_table.apply(
        lambda r: r["mirna_status"] if r["node_type"] == "miRNA" else "protein", axis=1
    )
    node_table.to_csv(NODES_CSV, index=False)
    print(node_table["color_group"].value_counts())


if __name__ == "__main__":
    main()
