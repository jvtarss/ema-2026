
import sqlite3, pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib.ticker as mticker

con = sqlite3.connect("ema_mirdeep2_union-definitivo-final.db")
df = pd.read_sql("""
SELECT mature_read_count, star_read_count FROM mirna_discovery_evidence
WHERE mature_read_count > 0 AND star_read_count > 0
""", con)

logx = np.log10(df.mature_read_count) - 1
logy = np.log10(df.star_read_count) - 1

slope, intercept = np.polyfit(logx, logy, 1)
xline = np.linspace(logx.min(), logx.max(), 100)
yline = slope * xline + intercept

def fmt(v, pos):
    return f"{v:.0f}"

fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(logx, logy, alpha=0.5)
ax.plot(xline, yline, color="red", linestyle="--", label=f"trend (slope={slope:.2f})")
ax.axhline(np.log10(10) - 1, color="gray", linestyle=":", linewidth=1)
ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt))
ax.set_xlabel("Mature read count (log10, 10¹ = 0)")
ax.set_ylabel("Star read count (log10, 10¹ = 0)")
ax.legend()
plt.tight_layout()
plt.savefig("figure3d_trend.png", dpi=300)
plt.savefig("figure3d_trend.pdf")
