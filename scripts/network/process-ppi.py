import pandas as pd

INPUT_TSV = "string_interactions_reciprocal.tsv"
OUTPUT_CSV = "ppi_edges_harmonized_dedup.csv"
CONFIDENCE_THRESHOLD = 0.400
OLD_PREFIX = "EUGRSUZ_"
NEW_PREFIX = "Eucgr."


def main():
    df = pd.read_csv(INPUT_TSV, sep="\t")
    df.columns = df.columns.str.lstrip("#")

    df_filtered = df[df["combined_score"] >= CONFIDENCE_THRESHOLD].copy()

    df_filtered["node1"] = df_filtered["node1"].str.replace(OLD_PREFIX, NEW_PREFIX, regex=False)
    df_filtered["node2"] = df_filtered["node2"].str.replace(OLD_PREFIX, NEW_PREFIX, regex=False)

    df_filtered["pair_key"] = df_filtered.apply(
        lambda r: tuple(sorted([r["node1"], r["node2"]])), axis=1
    )
    df_dedup = (
        df_filtered.sort_values("combined_score", ascending=False)
        .drop_duplicates(subset="pair_key", keep="first")
        .drop(columns="pair_key")
    )

    ppi_edges = df_dedup[["node1", "node2", "combined_score"]].rename(
        columns={"node1": "source", "node2": "target", "combined_score": "score"}
    )
    ppi_edges["interaction_type"] = "ppi-string"
    ppi_edges["inhibition_type"] = "Protein Interaction"
    ppi_edges = ppi_edges[["source", "target", "score", "interaction_type", "inhibition_type"]]

    ppi_edges.to_csv(OUTPUT_CSV, index=False)
    print(len(ppi_edges))


if __name__ == "__main__":
    main()
