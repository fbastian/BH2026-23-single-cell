"""Figure 2 - Human Liver Cell Atlas (HLiCA v1.0, CELLxGENE; Edgar, Portman et al. 2026) as used for SAMap (183 x 170 mm).
a deposited UMAP, 48 author_cell_type labels in the HLiCA authors' colours (colour_palette.R); b suspension; c study; d cells per study by suspension;
e HLiCA lineage composition per study (authors' lineage colours); f HLiCA canonical lineage markers (Zenodo 20380927) on the SAMap subsample (50,400 cells).
Outputs <D>/fig2_hlica.{pdf,png}."""
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
import numpy as np, pandas as pd
from matplotlib.lines import Line2D
from common import plt, MM, D, HLICA_COLORS, HLICA_LIN, HLICA_LIN_COLORS, HLICA_ORDER, ax_mm, letter, umap_ax
u = pd.read_parquet(f"{D}/data/hlica_umap.parquet"); U = u[["u1", "u2"]].values
u["lineage"] = u.label.map(HLICA_LIN)
col = HLICA_COLORS
fig = plt.figure(figsize=(183 * MM, 170 * MM))
ax = ax_mm(fig, 4, 5, 64, 64); umap_ax(ax, U, [col[l] for l in u.label], s=0.03)
ax.set_title(f"HLiCA v1.0: {len(u):,} cells and nuclei, {u.donor.nunique()} donors, {u.study.nunique()} studies", pad=2, loc="left")
for j, chunk in enumerate((HLICA_ORDER[:24], HLICA_ORDER[24:])):
    lg = ax_mm(fig, 70 + j * 31, 7, 30, 62); lg.axis("off")
    lg.legend([Line2D([], [], ls="", marker="o", ms=2.6, mfc=col[t], mec="none") for t in chunk], chunk, loc="upper left", frameon=False,
              handletextpad=0.1, labelspacing=0.28, borderaxespad=0, title="Cell type (author)" if j == 0 else " ", title_fontsize=5, alignment="left")
letter(fig, 1, 3, "a")
SUS = {"cell": "#56B4E9", "nucleus": "#D55E00"}
studies = u.study.value_counts().index.tolist()
STU = dict(zip(studies, ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]))
ax = ax_mm(fig, 136, 5, 26, 26); umap_ax(ax, U, [SUS[v] for v in u.suspension], s=0.02); ax.set_title("Suspension", pad=1)
lg = ax_mm(fig, 163, 12, 20, 10); lg.axis("off"); vc = u.suspension.value_counts()
lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=SUS[k], mec="none") for k in SUS], [f"{k}\n({vc[k]:,})" for k in SUS], frameon=False,
          loc="upper left", handletextpad=0.2, borderaxespad=0, labelspacing=0.4)
letter(fig, 133, 3, "b")
ax = ax_mm(fig, 136, 37, 26, 26); umap_ax(ax, U, [STU[v] for v in u.study], s=0.02); ax.set_title("Study", pad=1)
lg = ax_mm(fig, 136, 65, 47, 10); lg.axis("off")
lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=STU[k], mec="none") for k in studies], [k.replace("_", " ") for k in studies], frameon=False,
          loc="upper left", ncol=2, handletextpad=0.2, borderaxespad=0, columnspacing=0.6, labelspacing=0.25)
letter(fig, 133, 35, "c")
# d cells per study by suspension
ct = pd.crosstab(u.study, u.suspension).reindex(studies)
ax = ax_mm(fig, 30, 86, 50, 30); y = np.arange(len(studies))[::-1]
ax.barh(y, ct["cell"] / 1e3, color=SUS["cell"], height=0.7); ax.barh(y, ct["nucleus"] / 1e3, left=ct["cell"] / 1e3, color=SUS["nucleus"], height=0.7)
ax.set_yticks(y, [s.replace("_", " ") for s in studies]); ax.tick_params(axis="y", length=0, pad=1)
ax.set_xlabel("Cells and nuclei (thousands)")
dn = u.groupby("study").donor.nunique()
for yi, s in zip(y, studies): ax.text(ct.loc[s].sum() / 1e3 + 2, yi, f"{dn[s]} donors", va="center", fontsize=5)
ax.set_xlim(0, 175)
letter(fig, 1, 83, "d")
# e lineage composition per study
lc = pd.crosstab(u.study, u.lineage, normalize="index").reindex(studies)[[l for l in HLICA_LIN_COLORS if l in set(u.lineage)]]
ax = ax_mm(fig, 116, 86, 46, 30); left = np.zeros(len(studies))
for L in lc.columns:
    ax.barh(y, lc[L].values * 100, left=left, color=HLICA_LIN_COLORS[L], height=0.7, label=L); left += lc[L].values * 100
ax.set_yticks(y, [s.replace("_", " ") for s in studies]); ax.tick_params(axis="y", length=0, pad=1); ax.set_xlim(0, 100)
ax.set_xlabel("Share of study (%)")
lg = ax_mm(fig, 164, 86, 19, 30); lg.axis("off")
lg.legend([Line2D([], [], ls="", marker="s", ms=3, mfc=HLICA_LIN_COLORS[k], mec="none") for k in lc.columns], list(lc.columns), frameon=False, loc="upper left",
          handletextpad=0.2, borderaxespad=0, labelspacing=0.3, title="Lineage", title_fontsize=5, alignment="left")
letter(fig, 87, 83, "e")
# f lineage dot plot
mk = pd.read_csv(f"{D}/data/hlica_dot_genes.tsv", sep="\t"); mk = mk[mk["Marker Type"] == "Canonical"]
LORD = ["Hepatocytes", "Cholangiocyte", "Endothelia", "Mesenchyme", "Myeloid", "Lymphocytes"]
gl = []
for L in LORD: gl += mk[mk.Lineage == L].Gene.tolist()[:5] if L != "Lymphocytes" else [g for g in ["CD3E", "CD8A", "NKG7", "KLRF1", "MS4A1", "CD79A", "JCHAIN"] if g in mk.Gene.values]
mean = pd.read_csv(f"{D}/data/hlica_dot_mean.tsv", sep="\t", index_col=0); frac = pd.read_csv(f"{D}/data/hlica_dot_frac.tsv", sep="\t", index_col=0)
rows = [l for l in HLICA_LIN_COLORS if l in mean.index]
Mm = mean.loc[rows, gl]; Mm = (Mm - Mm.min()) / (Mm.max() - Mm.min()); F = frac.loc[rows, gl]
ax = ax_mm(fig, 26, 132, 120, 20)
xi, yi = np.meshgrid(np.arange(len(gl)), np.arange(len(rows))[::-1])
s_ = ax.scatter(xi.ravel(), yi.ravel(), s=F.values.ravel() * 9, c=Mm.values.ravel(), cmap="Reds", vmin=0, vmax=1, lw=0)
ax.set_xticks(range(len(gl)), gl, rotation=90, style="italic"); ax.set_yticks(range(len(rows))[::-1], rows)
ax.set_xlim(-0.7, len(gl) - 0.3); ax.set_ylim(-0.7, len(rows) - 0.3); ax.tick_params(length=0, pad=1)
for sp in ax.spines.values(): sp.set_visible(False)
xb = 0
for L in LORD:
    nL = sum(g in mk[mk.Lineage == L].Gene.values for g in gl)
    ax.text(xb + nL / 2 - 0.5, len(rows) - 0.1, {"Hepatocytes": "Hepatocyte", "Myeloid": "Myeloid Cell", "Lymphocytes": "Lymphocyte"}.get(L, L), ha="center", va="bottom", fontsize=5)
    ax.plot([xb - 0.4, xb + nL - 0.6], [len(rows) - 0.25] * 2, color="k", lw=0.5, clip_on=False); xb += nL
ax.set_title("HLiCA canonical lineage markers (SAMap subsample, 50,400 cells)", pad=9)
cax = ax_mm(fig, 152, 132, 2, 11); cb = fig.colorbar(s_, cax=cax); cb.set_label("Scaled mean", fontsize=5)
cb.ax.tick_params(labelsize=5, width=0.5, length=1.5); cb.outline.set_linewidth(0.5)
lg = ax_mm(fig, 164, 130, 18, 14); lg.axis("off")
lg.legend([plt.scatter([], [], s=f * 9, c="k", lw=0) for f in (0.25, 0.5, 1.0)], ["25", "50", "100"], title="% expressing", title_fontsize=5,
          frameon=False, loc="upper left", handletextpad=0.2, labelspacing=0.3, borderaxespad=0)
letter(fig, 1, 124, "f")
fig.savefig(f"{D}/fig2_hlica.pdf"); fig.savefig(f"{D}/fig2_hlica.png", dpi=300)
print("written")
