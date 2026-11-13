"""Figure 1 - Guilliams 2022 mouse liver atlas rebuilt from GEO (183 x 170 mm, Nature spec).
a UMAP annotation of record (TotalVI, strict QC); b suspension; c digestion; d cells kept per QC step (paper: 185,894);
e composition by type and suspension; f Table S1 top-3 markers per type (dot plot, mean log CP10k scaled per gene).
Inputs <D>/data (fig_00_data.py). Outputs <D>/fig1_mouse_atlas.{pdf,png}."""
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
import numpy as np, pandas as pd
from common import plt, MM, D, MOUSE_ORDER, MOUSE_LIN, shades, ax_mm, letter, umap_ax
from matplotlib.lines import Line2D
u = pd.read_parquet(f"{D}/data/mouse_aor_umap.parquet"); U = u[["u1", "u2"]].values
order = [t for t in MOUSE_ORDER if t in set(u.label)]
col = shades([t for t in order if t != "unresolved"], MOUSE_LIN); col["unresolved"] = (0.8, 0.8, 0.8)
fig = plt.figure(figsize=(183 * MM, 170 * MM))
# a
ax = ax_mm(fig, 4, 4, 68, 56); umap_ax(ax, U, [col[l] for l in u.label], s=0.08)
ax.set_title(f"Mouse liver, rebuilt: {len(u):,} cells and nuclei", pad=2)
n = u.label.value_counts()
lg = ax_mm(fig, 79, 4, 34, 56); lg.axis("off")
lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=col[t], mec="none") for t in order], [f"{t} ({n[t]:,})" for t in order],
          loc="upper left", frameon=False, handletextpad=0.2, labelspacing=0.25, borderaxespad=0, title="Cell type (n)", title_fontsize=5,
          alignment="left")
letter(fig, 1, 2, "a")
# b, c
SUS = {"cell": "#56B4E9", "nucleus": "#D55E00"}; DIG = {"Ex Vivo": "#009E73", "In Vivo": "#CC79A7", "Nuclei": "#000000"}
for (y, key, pal, ttl, L) in ((4, "suspension", SUS, "Suspension", "b"), (33, "digestion", DIG, "Digestion", "c")):
    ax = ax_mm(fig, 120, y, 27, 26); umap_ax(ax, U, [pal[v] for v in u[key]], s=0.04); ax.set_title(ttl, pad=1)
    lg = ax_mm(fig, 152, y + 6, 29, 16); lg.axis("off")
    vc = u[key].value_counts()
    lg.legend([Line2D([], [], ls="", marker="o", ms=3, mfc=pal[k], mec="none") for k in pal], [f"{k} ({vc.get(k, 0):,})" for k in pal],
              loc="upper left", frameon=False, handletextpad=0.2, borderaxespad=0)
    letter(fig, 116, y - 2, L)
# d QC flow
# counts read from the step outputs (review 2026-11-12: no typed-in numbers); strict removals in the order reasons were assigned
G = "/mnt/masld_vault/cross_species/guilliams2022"
st = pd.read_csv(f"{G}/01_raw/sample_table.tsv", sep="\t")
n_raw, n_qc = int(st.n_cells_file.sum()), int(pd.read_csv(f"{G}/06_doublets_ambient/per_sample_summary.tsv", sep="\t").n_kept.sum())
rr = pd.read_parquet(f"{G}/10_strict/cells_strict.parquet", columns=["strict_remove_reason"]).strict_remove_reason.astype(str).value_counts()
n = int(rr.sum()); steps = [(f"Fig 1B samples on GEO ({len(st)})", n_raw), ("Authors' MAD QC + Scrublet", n_qc), ("Doublet clusters removed", n)]
for r, lab in (("scDblFinder_doublet", "scDblFinder doublets"), ("dx_lt200", "< 200 UMI after decontX"), ("bridge", "Between-type (bridge) cells"),
               ("sort_mismatch", "Sort-gate mismatches (final)")):
    n -= int(rr.get(r, 0)); steps.append((lab, n))
assert n == int(rr["kept"]) == len(u), (n, rr["kept"], len(u))
ax = ax_mm(fig, 42, 67, 52, 38)
y = np.arange(len(steps))[::-1]
ax.barh(y, [v / 1e3 for _, v in steps], color=["#BBBBBB"] * 6 + ["#000000"], height=0.65)
for yi, (_, v) in zip(y, steps): ax.text(v / 1e3 + 4, yi, f"{v:,}", va="center", fontsize=5)
ax.axvline(185.894, ls="--", lw=0.5, color="#D55E00"); ax.text(185.894, len(steps) - 0.35, "paper 185,894", color="k", fontsize=5, ha="center")
ax.set_yticks(y, [s for s, _ in steps]); ax.set_xlabel("Cells and nuclei kept (thousands)"); ax.set_xlim(0, 470)
letter(fig, 1, 63, "d")
# e composition
comp = pd.crosstab(u.label, u.suspension).reindex(order)
ax = ax_mm(fig, 124, 67, 52, 38)
yy = np.arange(len(order))[::-1]
ax.barh(yy, comp["cell"] / 1e3, color=SUS["cell"], height=0.7, label="cell")
ax.barh(yy, comp["nucleus"] / 1e3, left=comp["cell"] / 1e3, color=SUS["nucleus"], height=0.7, label="nucleus")
ax.set_yticks(yy, order, fontsize=5); ax.tick_params(axis="y", length=0, pad=1); ax.set_xlabel("Cells and nuclei (thousands)")
ax.legend(frameon=False, loc="lower right", handlelength=1)
letter(fig, 101, 63, "e")
# f dot plot
genes = pd.read_csv(f"{D}/data/dotplot_genes.tsv", sep="\t", header=None, index_col=0)[1].str.split(",")
mean = pd.read_csv(f"{D}/data/dotplot_mean.tsv", sep="\t", index_col=0); frac = pd.read_csv(f"{D}/data/dotplot_frac.tsv", sep="\t", index_col=0)
rows = [t for t in order if t != "unresolved"]
gl = list(dict.fromkeys(g for t in rows for g in genes[t]))
M = mean.loc[rows, gl]; M = (M - M.min()) / (M.max() - M.min()); F = frac.loc[rows, gl]
ax = ax_mm(fig, 22, 117, 140, 33)
xi, yi = np.meshgrid(np.arange(len(gl)), np.arange(len(rows))[::-1])
sc_ = ax.scatter(xi.ravel(), yi.ravel(), s=(F.values.ravel() * 9), c=M.values.ravel(), cmap="Reds", vmin=0, vmax=1, lw=0)
ax.set_xticks(range(len(gl)), gl, rotation=90, fontsize=5, style="italic"); ax.set_yticks(range(len(rows))[::-1], rows, fontsize=5)
ax.set_xlim(-0.7, len(gl) - 0.3); ax.set_ylim(-0.7, len(rows) - 0.3); ax.tick_params(length=0, pad=1)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title("Table S1 top-3 marker genes per type (Guilliams et al. 2022)", pad=2)
cax = ax_mm(fig, 166, 119, 2, 12); cb = fig.colorbar(sc_, cax=cax); cb.set_label("Scaled mean", fontsize=5); cb.ax.tick_params(labelsize=5, width=0.5, length=1.5); cb.outline.set_linewidth(0.5)
lg = ax_mm(fig, 164, 135, 17, 14); lg.axis("off")
lg.legend([plt.scatter([], [], s=f * 9, c="k", lw=0) for f in (0.25, 0.5, 1.0)], ["25", "50", "100"], title="% expressing", title_fontsize=5,
          frameon=False, loc="upper left", handletextpad=0.2, labelspacing=0.3, borderaxespad=0)
letter(fig, 1, 112, "f")
fig.savefig(f"{D}/fig1_mouse_atlas.pdf"); fig.savefig(f"{D}/fig1_mouse_atlas.png", dpi=300)
print("written")
