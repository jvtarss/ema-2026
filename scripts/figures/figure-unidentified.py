import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

con = sqlite3.connect("ema_mirdeep2_union-definitivo-final.db")
df = pd.read_sql("""
SELECT mc.situation, mp.chr_scaf 
FROM mirna_core mc 
JOIN mirna_precursors mp ON mp.mirna_accession = mc.accession
""", con)
con.close()

df["chr_scaf"] = df["chr_scaf"].fillna("Unknown")
counts = df.groupby(["chr_scaf", "situation"]).size().unstack(fill_value=0)
counts = counts.sort_index()

if "known" not in counts.columns:
    counts["known"] = 0
if "novel" not in counts.columns:
    counts["novel"] = 0

fig, ax = plt.subplots(figsize=(14, 6))
x = range(len(counts.index))
width = 0.35

rects1 = ax.bar([i - width / 2 for i in x], counts["known"], width, label="known")
rects2 = ax.bar([i + width / 2 for i in x], counts["novel"], width, label="novel")

ax.set_xticks(x)
ax.set_xticklabels(counts.index, rotation=45, ha="right")
ax.set_ylabel("Number of miRNAs")
ax.legend()

for rect in rects1:
    height = rect.get_height()
    if height > 0:
        ax.annotate(f'{int(height)}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

for rect in rects2:
    height = rect.get_height()
    if height > 0:
        ax.annotate(f'{int(height)}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig("figure_regions_known_novel.png", dpi=300)
plt.savefig("figure_regions_known_novel.pdf")
