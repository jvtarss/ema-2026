
import sqlite3, itertools, pandas as pd, matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

con = sqlite3.connect("ema_mirdeep2_union-definitivo-final.db")
meta = pd.read_sql("""
SELECT sa.srr_accession, sa.tissue, sa.genotype, sa.condition, st.author_id
FROM samples sa JOIN studies st ON st.study_id = sa.study_id
""", con)
meta["group"] = (
    meta.author_id.fillna("") + "|" +
    meta.tissue.fillna("") + "|" +
    meta.genotype.fillna("") + "|" +
    meta.condition.fillna("")
)

expr = pd.read_sql("""
SELECT mirna_core_accession, srr_accession, cpm
FROM mirna_expression
WHERE cpm > 0
""", con)
expr = expr.merge(meta[["srr_accession", "group"]], on="srr_accession")

sets = {}
for g, sub in expr.groupby("group"):
    sets[g] = set(sub.mirna_core_accession)

names = list(sets.keys())
combos = []
for r in range(1, len(names) + 1):
    for combo in itertools.combinations(names, r):
        members = set.intersection(*[sets[n] for n in combo])
        excluded = set.union(*[sets[n] for n in names if n not in combo]) if len(combo) < len(names) else set()
        exclusive = members - excluded
        if len(exclusive) > 0:
            combos.append((combo, len(exclusive)))
combos.sort(key=lambda x: -x[1])
combos = combos[:30]

fig = plt.figure(figsize=(14, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0.05)
ax_bar = fig.add_subplot(gs[0])
ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

x = range(len(combos))
ax_bar.bar(x, [c[1] for c in combos], color="black")
for i, c in enumerate(combos):
    ax_bar.text(i, c[1] + 0.5, str(c[1]), ha="center", fontsize=7, rotation=90)
ax_bar.set_ylabel("Intersection size")
ax_bar.set_xticks([])

for i, (combo, _) in enumerate(combos):
    for j, name in enumerate(names):
        color = "black" if name in combo else "lightgray"
        ax_matrix.scatter(i, j, color=color, s=40)
for i, (combo, _) in enumerate(combos):
    idxs = [j for j, name in enumerate(names) if name in combo]
    if len(idxs) > 1:
        ax_matrix.plot([i, i], [min(idxs), max(idxs)], color="black", linewidth=1.5, zorder=0)
ax_matrix.set_yticks(range(len(names)))
ax_matrix.set_yticklabels(names, fontsize=7)
ax_matrix.set_xticks([])
ax_matrix.set_ylim(-0.5, len(names) - 0.5)

plt.tight_layout()
plt.savefig("figure_upset_all_samples_merged.png", dpi=300)
plt.savefig("figure_upset_all_samples_merged.pdf")
