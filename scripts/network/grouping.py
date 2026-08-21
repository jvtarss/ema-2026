import re
import pandas as pd

ENRICHMENT_FILE = "enrichment.all.tsv"
ENRICHMENT_SEP = "\t"

MASTER_EDGES_CSV = "master_edges.csv"

INTERACTIONS_FILE = "string_interactions_reciprocal.tsv"
INTERACTIONS_SEP = "\t"

PROTEIN_ALIASES_FILE = "string_protein_annotations.tsv"
PROTEIN_ALIASES_SEP = "\t"
ALIASES_ID_COLUMN = "node"
ALIASES_COLUMN = "other_names_and_aliases"

NODES_CSV = "node_attributes.csv"
NODES_ID_COLUMN = "node_id"

FDR_THRESHOLD = 0.05
GO_CATEGORIES_TO_USE = ["go process", "go function"]
FALLBACK_LABEL = "Other"

OLD_PREFIX = "EUGRSUZ_"
NEW_PREFIX = "Eucgr."
TAXON_PREFIX = "71139."

EXCLUDED_TERM_SUBSTRINGS = [
    "cellular anatomical entity",
    "organelle",
    "biological_process",
    "metabolic process",
    "binding",
]

CATEGORY_COL_CANDIDATES = ["category"]
TERM_COL_CANDIDATES = ["term description", "term_description", "description"]
FDR_COL_CANDIDATES = ["false discovery rate", "fdr", "false_discovery_rate"]
STRENGTH_COL_CANDIDATES = ["strength"]
GENES_COL_CANDIDATES = [
    "matching proteins in your network (ids)",
    "matching_proteins_in_your_network_ids",
    "matching proteins in your network (labels)",
]

INTERACTIONS_NODE_COLS = [("node1", "node1_string_id"), ("node2", "node2_string_id")]


def find_column(df, candidates, label, required=True):
    for c in candidates:
        for col in df.columns:
            if col.strip().lower() == c.strip().lower():
                return col
    if required:
        raise KeyError(f"Nenhuma coluna encontrada para '{label}'. Colunas disponiveis: {list(df.columns)}")
    return None


def is_excluded_term(description):
    lowered = description.lower()
    return any(sub in lowered for sub in EXCLUDED_TERM_SUBSTRINGS)


def strip_taxon_prefix(protein_id):
    if protein_id.startswith(TAXON_PREFIX):
        return protein_id[len(TAXON_PREFIX):]
    return protein_id


def add_to_map(mapping, key, value):
    mapping.setdefault(key, set()).add(value)


def build_eucgr_map_from_interactions(path):
    mapping = {}
    df = pd.read_csv(path, sep=INTERACTIONS_SEP)
    df.columns = df.columns.str.lstrip("#").str.strip()

    for name_col, id_col in INTERACTIONS_NODE_COLS:
        if name_col not in df.columns or id_col not in df.columns:
            continue
        for _, row in df.iterrows():
            name_value = str(row[name_col])
            id_value = str(row[id_col])
            if name_value.startswith(OLD_PREFIX):
                eucgr_id = NEW_PREFIX + name_value[len(OLD_PREFIX):]
                raw_id = strip_taxon_prefix(id_value)
                add_to_map(mapping, raw_id, eucgr_id)
                add_to_map(mapping, eucgr_id, raw_id)
    return mapping


def build_eucgr_map_from_aliases(path, mapping):
    ann = pd.read_csv(path, sep=PROTEIN_ALIASES_SEP)
    ann.columns = ann.columns.str.lstrip("#")
    id_col = find_column(ann, [ALIASES_ID_COLUMN], "node id", required=False)
    alias_col = find_column(ann, [ALIASES_COLUMN], "aliases", required=False)
    if id_col is None or alias_col is None:
        return mapping

    for _, row in ann.iterrows():
        raw_node = row[id_col]
        aliases_text = str(row[alias_col])
        for token in aliases_text.split(","):
            token = token.strip()
            if token.startswith(OLD_PREFIX):
                eucgr_id = NEW_PREFIX + token[len(OLD_PREFIX):]
                add_to_map(mapping, raw_node, eucgr_id)
    return mapping


def harmonize_all(protein_id, eucgr_map):
    ids = {protein_id}
    if protein_id in eucgr_map:
        ids |= eucgr_map[protein_id]
    if OLD_PREFIX in protein_id:
        ids.add(protein_id.replace(OLD_PREFIX, NEW_PREFIX))
    return ids


def load_go_terms(eucgr_map):
    df = pd.read_csv(ENRICHMENT_FILE, sep=ENRICHMENT_SEP)
    df.columns = df.columns.str.lstrip("#").str.strip()

    category_col = find_column(df, CATEGORY_COL_CANDIDATES, "category")
    term_col = find_column(df, TERM_COL_CANDIDATES, "term description")
    fdr_col = find_column(df, FDR_COL_CANDIDATES, "FDR")
    strength_col = find_column(df, STRENGTH_COL_CANDIDATES, "strength")
    genes_col = find_column(df, GENES_COL_CANDIDATES, "matching proteins")

    df = df[df[category_col].str.strip().str.lower().isin(GO_CATEGORIES_TO_USE)]

    df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")
    df[strength_col] = pd.to_numeric(df[strength_col], errors="coerce")
    df = df.dropna(subset=[fdr_col, strength_col])
    df = df[df[fdr_col] <= FDR_THRESHOLD]
    df = df[~df[term_col].apply(is_excluded_term)]
    df = df.sort_values([strength_col, fdr_col], ascending=[False, True])

    terms = []
    for _, row in df.iterrows():
        raw_proteins = [g.strip() for g in str(row[genes_col]).split(",") if g.strip()]
        proteins = set()
        for p in raw_proteins:
            p = strip_taxon_prefix(p)
            proteins |= harmonize_all(p, eucgr_map)
        terms.append((str(row[term_col]), row[fdr_col], proteins))
    return terms


def load_mirna_targets():
    df = pd.read_csv(MASTER_EDGES_CSV)
    df = df[df["interaction_type"] == "miRNA-target"]

    mirna_targets = {}
    for _, row in df.iterrows():
        mirna_targets.setdefault(row["source"], set()).add(row["target"])
    return mirna_targets


def sanitize_term(text):
    return text.replace(",", ";")


def assign_protein_groups(go_terms):
    protein_group = {}
    protein_group_fdr = {}

    for term_description, fdr_value, term_proteins in go_terms:
        for protein in term_proteins:
            if protein in protein_group:
                continue
            protein_group[protein] = sanitize_term(term_description)
            protein_group_fdr[protein] = fdr_value

    return protein_group, protein_group_fdr


def main():
    eucgr_map = build_eucgr_map_from_interactions(INTERACTIONS_FILE)
    print(f"{len(eucgr_map)} chaves no mapa vindas de string_interactions_reciprocal.tsv")

    eucgr_map = build_eucgr_map_from_aliases(PROTEIN_ALIASES_FILE, eucgr_map)
    print(f"{len(eucgr_map)} chaves no mapa apos somar aliases de string_protein_annotations.tsv")

    go_terms = load_go_terms(eucgr_map)
    print(f"{len(go_terms)} termos GO Process/Function significativos apos filtros")

    mirna_targets = load_mirna_targets()
    print(f"{len(mirna_targets)} miRNAs com alvo mapeado em master_edges.csv")

    total_protein_universe = set()
    for _, _, proteins in go_terms:
        total_protein_universe |= proteins
    total_target_universe = set()
    for targets in mirna_targets.values():
        total_target_universe |= targets
    overlap = total_protein_universe & total_target_universe
    print(f"sobreposicao bruta entre universo de proteinas do enriquecimento e universo de alvos de miRNA: {len(overlap)}")

    mirna_group = {}
    mirna_group_fdr = {}
    mirna_group_hits = {}

    for mirna_id, targets in mirna_targets.items():
        assigned = False
        for term_description, fdr_value, term_proteins in go_terms:
            hits = targets & term_proteins
            if hits:
                mirna_group[mirna_id] = sanitize_term(term_description)
                mirna_group_fdr[mirna_id] = fdr_value
                mirna_group_hits[mirna_id] = len(hits)
                assigned = True
                break
        if not assigned:
            mirna_group[mirna_id] = FALLBACK_LABEL
            mirna_group_fdr[mirna_id] = None
            mirna_group_hits[mirna_id] = 0

    protein_group, protein_group_fdr = assign_protein_groups(go_terms)

    nodes = pd.read_csv(NODES_CSV)
    nodes["target_go_term"] = nodes[NODES_ID_COLUMN].map(mirna_group)
    nodes["target_go_term_fdr"] = nodes[NODES_ID_COLUMN].map(mirna_group_fdr)
    nodes["target_go_term_hits"] = nodes[NODES_ID_COLUMN].map(mirna_group_hits)
    nodes["protein_go_term"] = nodes[NODES_ID_COLUMN].map(protein_group)
    nodes["protein_go_term_fdr"] = nodes[NODES_ID_COLUMN].map(protein_group_fdr)

    nodes["color_group"] = nodes.apply(
        lambda r: r["mirna_status"] if r["node_type"] == "miRNA"
        else protein_group.get(r[NODES_ID_COLUMN], FALLBACK_LABEL),
        axis=1,
    )

    text_columns = nodes.select_dtypes(include="object").columns
    for col in text_columns:
        if col == NODES_ID_COLUMN:
            continue
        nodes[col] = nodes[col].astype(str).str.replace(",", ";", regex=False)
        nodes[col] = nodes[col].replace("nan", "")

    nodes.to_csv(NODES_CSV, index=False)

    print("color_group (miRNA known/novel + protein GO group, coluna unica):")
    print(nodes["color_group"].value_counts())
    print("\ncolor_group restrito a miRNA:")
    print(nodes.loc[nodes["node_type"] == "miRNA", "color_group"].value_counts())
    print("\ncolor_group restrito a protein:")
    print(nodes.loc[nodes["node_type"] == "protein", "color_group"].value_counts())


if __name__ == "__main__":
    main()
