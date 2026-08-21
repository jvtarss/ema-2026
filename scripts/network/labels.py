import re
import pandas as pd
from collections import Counter

ANNOTATIONS_FILE = "string_protein_annotations.tsv"
ANNOTATIONS_SEP = "\t"
ID_COLUMN = "node"
ALIASES_COLUMN = "other_names_and_aliases"

NODES_CSV = "node_attributes.csv"
NODES_ID_COLUMN = "node_id"

EUGRSUZ_PATTERN = re.compile(r"EUGRSUZ_([A-Z]\d{5})")
NEW_PREFIX = "Eucgr."

ID_LIKE_PATTERNS = [
    re.compile(r"^\d+$"),
    re.compile(r"^EUGRSUZ_"),
    re.compile(r"^A0A[0-9A-Z]+"),
    re.compile(r"^[A-Z]{1,3}_\d+"),
    re.compile(r"^egr:"),
    re.compile(r"^[A-Z0-9]+_EUCGR$"),
]

WHITELIST = {
    "nac domain": "NAC",
    "nb-arc": "NBARC",
    "tir domain": "TIR",
    "leucine-rich repeat": "LRR",
    "sbp domain": "SBP",
    "grf": "GRF",
    "myb": "MYB",
    "wrky": "WRKY",
    "bhlh": "bHLH",
    "bzip": "bZIP",
    "ank_rep_region": "ANK",
    "wd_repeats_region": "WD",
    "f-box": "FBOX",
    "zinc finger": "ZNF",
    "ring finger": "RING",
    "heat shock": "HSP",
    "cytochrome p450": "CYP450",
    "gtpase": "GTPase",
    "protein kinase": "KIN",
    "lrrnt": "LRRNT",
    "c2 domain": "C2",
    "cid domain": "CID",
    "dna primase": "PRIM",
}

STOPWORDS = {
    "a", "an", "the", "of", "and", "or", "to", "in", "type", "containing",
    "like", "related", "associated", "putative", "probable", "predicted",
    "uncharacterized",
}


def detect_description_column(df):
    candidates = [c for c in df.columns if c not in (ID_COLUMN, "identifier", ALIASES_COLUMN)]
    avg_word_counts = {
        c: df[c].astype(str).str.split().str.len().mean() for c in candidates
    }
    return max(avg_word_counts, key=avg_word_counts.get)


def clean_description(text):
    if not isinstance(text, str) or not text.strip():
        return ""
    text = text.split(";")[0]
    text = text.split(". ")[0]
    text = text.strip().rstrip(".")
    return text


def match_whitelist(text):
    lowered = text.lower()
    for key, label in WHITELIST.items():
        if key in lowered:
            return label
    return None


def build_fallback_label(text):
    parts = []
    for word in text.split():
        subtokens = word.split("-")
        for sub in subtokens:
            sub_clean = sub.strip()
            if not sub_clean:
                continue
            if sub_clean.lower() in STOPWORDS:
                continue
            if len(sub_clean) <= 3 and any(c.isupper() for c in sub_clean):
                parts.append(sub_clean)
            else:
                parts.append(sub_clean[0].upper())
    return "".join(parts) if parts else ""


def build_label(fallback_id, description):
    cleaned = clean_description(description)
    if not cleaned:
        return fallback_id

    whitelist_hit = match_whitelist(cleaned)
    if whitelist_hit:
        return whitelist_hit

    fallback = build_fallback_label(cleaned)
    return fallback if fallback else fallback_id


def extract_eucgr_alias(aliases_text):
    if not isinstance(aliases_text, str):
        return None
    match = EUGRSUZ_PATTERN.search(aliases_text)
    if match:
        return NEW_PREFIX + match.group(1)
    return None


def is_id_like(token):
    if " " not in token:
        for pattern in ID_LIKE_PATTERNS:
            if pattern.match(token):
                return True
        if len(token) < 4:
            return True
    return False


def extract_common_name(aliases_text, fallback_description):
    if not isinstance(aliases_text, str):
        return fallback_description

    tokens = [t.strip() for t in aliases_text.split(",") if t.strip()]
    candidates = [t for t in tokens if " " in t and not is_id_like(t)]

    if not candidates:
        return fallback_description

    return max(candidates, key=len)


def main():
    ann = pd.read_csv(ANNOTATIONS_FILE, sep=ANNOTATIONS_SEP)
    ann.columns = ann.columns.str.lstrip("#")

    description_col = detect_description_column(ann)
    print("coluna de descricao detectada:", description_col)

    label_map = {}
    description_map = {}
    common_name_map = {}
    for _, row in ann.iterrows():
        raw_node = row[ID_COLUMN]
        description = row[description_col]
        cleaned = clean_description(description)
        label = build_label(raw_node, description)
        common_name = extract_common_name(row[ALIASES_COLUMN], cleaned)

        label_map[raw_node] = label
        description_map[raw_node] = cleaned
        common_name_map[raw_node] = common_name

        eucgr_id = extract_eucgr_alias(row[ALIASES_COLUMN])
        if eucgr_id:
            label_map[eucgr_id] = label
            description_map[eucgr_id] = cleaned
            common_name_map[eucgr_id] = common_name

    nodes = pd.read_csv(NODES_CSV)
    raw_labels = nodes[NODES_ID_COLUMN].map(label_map)
    raw_labels = raw_labels.fillna(nodes[NODES_ID_COLUMN])

    used = Counter()
    final = []
    for lbl in raw_labels:
        used[lbl] += 1
        final.append(lbl if used[lbl] == 1 else f"{lbl}-{used[lbl]}")

    nodes["display_label"] = final
    nodes["string_annotation"] = nodes[NODES_ID_COLUMN].map(description_map).fillna("")
    nodes["common_name"] = nodes[NODES_ID_COLUMN].map(common_name_map).fillna("")
    nodes.to_csv(NODES_CSV, index=False)

    n_found = (nodes[NODES_ID_COLUMN] != raw_labels).sum()
    n_missing = (nodes[NODES_ID_COLUMN] == raw_labels).sum()
    print(nodes["display_label"].nunique(), "labels unicos gerados")
    print(n_found, "nos com anotacao encontrada")
    print(n_missing, "nos sem anotacao encontrada")


if __name__ == "__main__":
    main()
