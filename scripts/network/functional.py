import re
import pandas as pd
from collections import Counter

ENRICHMENT_FILE = "enrichment.all.tsv"
ENRICHMENT_SEP = "\t"

PROTEIN_ALIASES_FILE = "string_protein_annotations.tsv"
PROTEIN_ALIASES_SEP = "\t"
ALIASES_ID_COLUMN = "node"
ALIASES_COLUMN = "other_names_and_aliases"

NODES_CSV = "node_attributes.csv"
NODES_ID_COLUMN = "node_id"

FDR_THRESHOLD = 0.05
TOP_UNCLASSIFIED_TERMS_TO_SHOW = 25

OLD_PREFIX = "EUGRSUZ_"
NEW_PREFIX = "Eucgr."
EUGRSUZ_PATTERN = re.compile(r"EUGRSUZ_([A-Z]\d{5})")
TAXON_PREFIX = "71139."

CATEGORY_COL_CANDIDATES = ["category"]
TERM_COL_CANDIDATES = ["term description", "term_description", "description"]
FDR_COL_CANDIDATES = ["false discovery rate", "fdr", "false_discovery_rate"]
GENES_COL_CANDIDATES = [
    "matching proteins in your network (ids)",
    "matching_proteins_in_your_network_ids",
    "matching proteins in your network (labels)",
]

ONTOLOGY_PRIORITY = [
    "pfam",
    "interpro",
    "smart",
    "rctm",
    "keyword",
    "networkneighboral",
    "go function",
    "go process",
    "go component",
    "compartments",
]

EXCLUDED_TERM_SUBSTRINGS = [
    "cellular anatomical entity",
    "organelle",
    "cellular process",
    "biological_process",
    "metabolic process",
]

CATEGORY_KEYWORDS = [
    ("Plant immune receptor signaling (NBS-LRR / TIR)", [
        "disease resistance", "resistance protein", "nbs-lrr", "nb-arc",
        "tir domain", "toll", "interleukin", "leucine-rich repeat",
        "defense response", "hypersensitive response", "pathogen",
        "immune", "innate immunity", "r protein", "avirulence", "neutrophil",
    ]),
    ("Cell wall, lignin & secondary metabolism", [
        "cell wall", "lignin", "cellulose", "hemicellulose", "xylan",
        "phenylpropanoid", "flavonoid", "pectin", "glycosyltransferase",
        "laccase", "peroxidase", "cinnamyl alcohol dehydrogenase",
        "caffeoyl", "monolignol", "terpenoid", "wax biosynthesis", "cutin",
        "secondary metabolite",
    ]),
    ("Hormone signaling & growth control", [
        "auxin", "hormone", "signal transduction", "kinase signaling",
        "gibberellin", "cytokinin", "ethylene", "brassinosteroid",
        "jasmonate", "abscisic acid", "growth regulation",
        "auxin response factor", "growth factor", "receptor kinase",
        "g-protein",
    ]),
    ("Transcriptional & developmental regulators", [
        "transcription factor", "transcription regulation", "dna-binding",
        "zinc finger", "homeobox", "myb", "wrky", "bhlh", "bzip",
        "nac domain", "sbp domain", "squamosa promoter",
        "developmental process", "meristem", "phase change", "embryogenesis",
        "flowering", "floral", "seed development", "root development",
        "shoot development",
    ]),
    ("Photosynthesis, redox & stress response", [
        "photosynthesis", "chlorophyll", "chloroplast", "thylakoid",
        "redox", "oxidative stress", "reactive oxygen", "antioxidant",
        "catalase", "superoxide dismutase", "glutathione", "abiotic stress",
        "cold stress", "drought", "heat stress",
    ]),
    ("DNA replication, repair & chromatin", [
        "dna replication", "dna repair", "chromatin", "histone",
        "nucleosome", "helicase", "primase", "dna polymerase",
        "cell cycle", "mitosis", "meiosis",
    ]),
    ("Transport & membrane trafficking", [
        "transporter activity", "vesicle", "endocytosis", "exocytosis",
        "membrane trafficking", "golgi", "endoplasmic reticulum",
        "ion transport", "abc transporter", "aquaporin", "channel activity",
        "transport",
    ]),
    ("RNA processing, translation & protein turnover", [
        "translation", "ribosom", "trna", "rna binding", "methyltransferase",
        "splicing", "proteolysis", "ubiquitin", "protein modification",
        "rna processing", "spliceosome", "proteasome", "chaperone",
        "heat shock protein", "protein folding",
    ]),
    ("Primary metabolism & biosynthesis", [
        "biosynthetic process", "catabolic process", "oxidoreductase",
        "transferase activity", "hydrolase activity", "amino acid metabolism",
        "lipid metabolism", "carbohydrate metabolism", "nucleotide metabolism",
        "fatty acid", "sugar metabolism", "starch", "biosynthesis",
        "metabolism", "enzyme activity",
    ]),
]

FALLBACK_CATEGORY = "Other / unclassified"


def find_column(df, candidates, label, required=True):
    for c in candidates:
        for col in df.columns:
            if col.strip().lower() == c.strip().lower():
                return col
    if required:
        raise KeyError(f"Nenhuma coluna encontrada para '{label}'. Colunas disponiveis: {list(df.columns)}")
    return None


def ontology_rank(category_value):
    lowered = category_value.strip().lower()
    for i, tag in enumerate(ONTOLOGY_PRIORITY):
        if tag in lowered:
            return i
    return len(ONTOLOGY_PRIORITY)


def is_excluded_term(description):
    lowered = description.lower()
    return any(sub in lowered for sub in EXCLUDED_TERM_SUBSTRINGS)


def classify_term(description):
    lowered = description.lower()
    for category, keywords in CATEGORY_KEYWORDS:
        for kw in keywords:
            if kw in lowered:
                return category
    return FALLBACK_CATEGORY


def build_eucgr_map(path):
    ann = pd.read_csv(path, sep=PROTEIN_ALIASES_SEP)
    ann.columns = ann.columns.str.lstrip("#")
    id_col = find_column(ann, [ALIASES_ID_COLUMN], "node id", required=False)
    alias_col = find_column(ann, [ALIASES_COLUMN], "aliases", required=False)
    if id_col is None or alias_col is None:
        return {}

    mapping = {}
    for _, row in ann.iterrows():
        raw_node = row[id_col]
        match = EUGRSUZ_PATTERN.search(str(row[alias_col]))
        if match:
            mapping[raw_node] = NEW_PREFIX + match.group(1)
    return mapping


def strip_taxon_prefix(protein_id):
    if protein_id.startswith(TAXON_PREFIX):
        return protein_id[len(TAXON_PREFIX):]
    return protein_id


def harmonize(protein_id, eucgr_map):
    if protein_id in eucgr_map:
        return eucgr_map[protein_id]
    if OLD_PREFIX in protein_id:
        return protein_id.replace(OLD_PREFIX, NEW_PREFIX)
    return protein_id


def assign_categories(term_records, eucgr_map):
    protein_category = {}
    protein_term = {}
    protein_fdr = {}
    protein_all_terms = {}

    for term_description, fdr_value, proteins in term_records:
        broad_category = classify_term(term_description)
        for protein in proteins:
            harmonized = harmonize(protein, eucgr_map)
            candidates = {protein, harmonized}

            for c in candidates:
                protein_all_terms.setdefault(c, []).append(term_description)

            if any(c in protein_category for c in candidates):
                continue
            for c in candidates:
                protein_category[c] = broad_category
                protein_term[c] = term_description
                protein_fdr[c] = fdr_value

    return protein_category, protein_term, protein_fdr, protein_all_terms


def main():
    df = pd.read_csv(ENRICHMENT_FILE, sep=ENRICHMENT_SEP)
    df.columns = df.columns.str.lstrip("#").str.strip()

    category_col = find_column(df, CATEGORY_COL_CANDIDATES, "category")
    term_col = find_column(df, TERM_COL_CANDIDATES, "term description")
    fdr_col = find_column(df, FDR_COL_CANDIDATES, "FDR")
    genes_col = find_column(df, GENES_COL_CANDIDATES, "matching proteins")

    df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")
    df = df.dropna(subset=[fdr_col])
    df = df[df[fdr_col] <= FDR_THRESHOLD]

    df = df[~df[term_col].apply(is_excluded_term)]

    df["ontology_rank"] = df[category_col].apply(ontology_rank)
    df = df.sort_values(["ontology_rank", fdr_col], ascending=[True, True])

    eucgr_map = build_eucgr_map(PROTEIN_ALIASES_FILE)

    term_records = []
    for _, row in df.iterrows():
        raw_proteins = [g.strip() for g in str(row[genes_col]).split(",") if g.strip()]
        proteins = [strip_taxon_prefix(p) for p in raw_proteins]
        term_records.append((str(row[term_col]), row[fdr_col], proteins))

    protein_category, protein_term, protein_fdr, protein_all_terms = assign_categories(term_records, eucgr_map)

    nodes = pd.read_csv(NODES_CSV)
    nodes["functional_category"] = nodes[NODES_ID_COLUMN].map(protein_category).fillna(FALLBACK_CATEGORY)
    nodes["functional_term_description"] = nodes[NODES_ID_COLUMN].map(protein_term).fillna("")
    nodes["functional_term_fdr"] = nodes[NODES_ID_COLUMN].map(protein_fdr)
    nodes.to_csv(NODES_CSV, index=False)

    print(nodes["functional_category"].value_counts())

    unclassified_ids = nodes[nodes["functional_category"] == FALLBACK_CATEGORY][NODES_ID_COLUMN]
    leftover_terms = Counter()
    for node_id in unclassified_ids:
        for term in protein_all_terms.get(node_id, []):
            leftover_terms[term] += 1

    if leftover_terms:
        print(f"\nTermos mais frequentes entre os {FALLBACK_CATEGORY} (nao capturados por nenhuma keyword):")
        for term, count in leftover_terms.most_common(TOP_UNCLASSIFIED_TERMS_TO_SHOW):
            print(f"  {count:4d}  {term}")


if __name__ == "__main__":
    main()
