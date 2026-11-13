"""Figure 3 - SAMap human (HLiCA v1.0) x mouse (Guilliams 2022, rebuilt) liver (183 x 170 mm, Nature spec).
a,b unintegrated baseline (same cells, MGI 1:1 orthologs, normalise + log1p, 3,000 HVG, PCA 50, UMAP; fig_00_data.py) by species / lineage;
c,d SAMap joint UMAP (sc-samap 3.0.1) by species / lineage; e best matches mouse type -> HLiCA type (mapping score >= 0.4, width = score);
f mapping scores for all mouse x HLiCA types (lineage bars on the margins).
Note: HLiCA contains the Guilliams 2022 human samples (51,554 cells); the mouse atlas is the same lab's mouse data.
Layout rule (2026-11-12 user): no overlapping text - every label block has its own reserved mm band.
Gene pairs behind each match: separate Table 1 (tab1_data.py + tab1_kable.R; user 2026-11-12, panel removed, a-d enlarged).
Outputs <D>/fig3_samap.{pdf,png}."""
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
import numpy as np, pandas as pd
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D
from common import plt, MM, D, WONG, LIN_ORDER, MOUSE_LIN, HUMAN_LIN, MOUSE_ORDER, HUMAN_ORDER, ax_mm, letter, umap_ax
O = "/mnt/masld_vault/cross_species/samap_run/out"
SPC = {"human": "#D55E00", "mouse": "#0072B2"}
def lin_of(sp, l): return (HUMAN_LIN if sp == "human" else MOUSE_LIN).get(l, "Other")
fig = plt.figure(figsize=(183 * MM, 170 * MM))
# a-d: 2 x 2 UMAPs, unintegrated row first (y 5-38), SAMap row (y 45-78)
for i, (f, ttl) in enumerate((("unintegrated_umap", "Unintegrated"), ("samap_umap", "SAMap"))):
    u = pd.read_parquet(f"{D}/data/{f}.parquet"); U = u[["u1", "u2"]].values
    y = 7 + i * 49
    ax = ax_mm(fig, 5, y, 42, 42); umap_ax(ax, U, [SPC[s] for s in u.species], s=0.05); ax.set_title(f"{ttl}: species", pad=1)
    ax = ax_mm(fig, 50, y, 42, 42); umap_ax(ax, U, [WONG[lin_of(s, l)] for s, l in zip(u.species, u.label)], s=0.05)
    ax.set_title(f"{ttl}: lineage", pad=1)
    letter(fig, 1, y - 4, "ac"[i]); letter(fig, 46, y - 4, "bd"[i])
n = pd.read_parquet(f"{D}/data/samap_umap.parquet").species.value_counts()
lg = ax_mm(fig, 93, 10, 14, 10); lg.axis("off")
lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=SPC[k], mec="none") for k in SPC], [f"{k}\n({n[k]:,})" for k in SPC], title="Species",
          title_fontsize=5, frameon=False, loc="upper left", handletextpad=0.2, borderaxespad=0, labelspacing=0.5, alignment="left")
lg = ax_mm(fig, 93, 34, 14, 30); lg.axis("off")
lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=WONG[k], mec="none") for k in LIN_ORDER], LIN_ORDER, title="Lineage", title_fontsize=5,
          frameon=False, loc="upper left", handletextpad=0.2, borderaxespad=0, labelspacing=0.3, alignment="left")
# mapping table
M = pd.read_csv(f"{O}/mapping_table.tsv", sep="\t", index_col=0)
mrows = [t for t in MOUSE_ORDER if "mm_" + t in M.index]; hcols = [t for t in HUMAN_ORDER if "hs_" + t in M.columns]
T = pd.DataFrame(M.loc[["mm_" + t for t in mrows], ["hs_" + t for t in hcols]].values, index=mrows, columns=hcols)
# e: x bands - left labels 104-123, curves 124-142, right labels 143-183; y 6-80
thr = 0.4
E = [(m, h, T.loc[m, h]) for m in mrows for h in hcols if T.loc[m, h] >= thr]
hs_used = [h for h in hcols if any(e[1] == h for e in E)]
ax = ax_mm(fig, 126, 8, 16, 90); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
yl = {m: 1 - (i + 0.5) / len(mrows) for i, m in enumerate(mrows)}; yr = {h: 1 - (i + 0.5) / len(hs_used) for i, h in enumerate(hs_used)}
for m, h, s in sorted(E, key=lambda e: e[2]):
    p = Path([(0, yl[m]), (0.5, yl[m]), (0.5, yr[h]), (1, yr[h])], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(p, fc="none", ec=WONG[MOUSE_LIN[m]], lw=0.3 + 1.8 * (s - thr) / (1 - thr), alpha=0.8))
for m in mrows: ax.text(-0.05, yl[m], m, ha="right", va="center", fontsize=5)
for h in hs_used: ax.text(1.05, yr[h], h, ha="left", va="center", fontsize=5)
fig.text(125 / 183, 1 - 5.5 / 170, "Mouse", ha="right", va="bottom", fontsize=6)
fig.text(143 / 183, 1 - 5.5 / 170, "Human (HLiCA)", ha="left", va="bottom", fontsize=6)
fig.text(134 / 183, 1 - 5.5 / 170, f"score ≥ {thr}", ha="center", va="bottom", fontsize=5)
letter(fig, 105, 3, "e")
# g heatmap: y 107-139, column labels below to ~166
ax = ax_mm(fig, 24, 107, 133, 37)
im = ax.imshow(T.values, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
ax.set_yticks(range(len(mrows)), mrows, fontsize=5); ax.set_xticks(range(len(hcols)), hcols, rotation=90, fontsize=5)
ax.tick_params(length=0, pad=1)
for sp in ax.spines.values(): sp.set_visible(False)
xb = 0
for L in LIN_ORDER:
    nL = sum(HUMAN_LIN[h] == L for h in hcols)
    if nL: ax.plot([xb - 0.45, xb + nL - 0.55], [-0.95, -0.95], color=WONG[L], lw=1.6, clip_on=False, solid_capstyle="butt"); xb += nL
yb = 0
for L in LIN_ORDER:
    nL = sum(MOUSE_LIN[m] == L for m in mrows)
    if nL: ax.plot([len(hcols) - 0.5 + 0.7] * 2, [yb - 0.45, yb + nL - 0.55], color=WONG[L], lw=1.6, clip_on=False, solid_capstyle="butt"); yb += nL
ax.set_ylabel("Mouse", labelpad=1); ax.set_title("SAMap mapping scores, HLiCA types (columns)", pad=5)
cax = ax_mm(fig, 165, 110, 2, 14); cb = fig.colorbar(im, cax=cax); cb.set_label("Mapping score", fontsize=5)
cb.ax.tick_params(labelsize=5, width=0.5, length=1.5); cb.outline.set_linewidth(0.5)
letter(fig, 1, 104, "f")
fig.savefig(f"{D}/fig3_samap.pdf"); fig.savefig(f"{D}/fig3_samap.png", dpi=300)
print("edges", len(E), "human types in e", len(hs_used))
