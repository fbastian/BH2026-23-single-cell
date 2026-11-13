"""Step 06c - assemble scDblFinder + decontX results for the draft-QC cells (328,041).

Doublets (MASLD rule): consensus_doublet = Scrublet AND scDblFinder; removal only on consensus. Per-cluster (TotalVI Leiden r3)
doublet rates are reported for the authors' style cluster-level removal (decided by hand, not here).
Ambient: decontX counts per sample; Alb diagnostic = Alb CP10k in non-epithelial draft types before/after, Alb kept in
hepatocytes.
Outputs <OUT>/06_doublets_ambient/:
  guilliams2022_fig1b_decontx.h5ad  X = decontX counts, layers['counts'] = raw counts, obs flags, obsm adt / adt_measured
  cells_doublets_ambient.parquet, per_sample_summary.tsv, per_cluster_r3_summary.tsv, figures/*.png|pdf
"""
import glob
import os
import sys

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/06_doublets_ambient"
os.makedirs(f"{D}/figures", exist_ok=True)

a = ad.read_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad")
f = pd.read_parquet(f"{OUT}/02_qc/guilliams2022_fig1b_qc_flags.parquet")
E = pd.read_parquet(f"{OUT}/05_annotation_draft/cells_annotation_draft.parquet")
assert (E.index == a.obs_names).all()
res = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{D}/per_sample/*.cells.parquet")])
assert len(res) == 389056 and res.index.is_unique

# decontX counts, kept cells in each sample's order (= raw object order = 02_qc order within sample)
blocks, names = [], []
for s, idx in a.obs.groupby("sample", observed=True).indices.items():
    m = sp.load_npz(f"{D}/per_sample/{s}.decontx.npz")
    assert m.shape[0] == len(idx), s
    blocks.append(m)
    names.append(a.obs_names[idx])
order = np.concatenate([a.obs_names.get_indexer(n) for n in names])
M = sp.vstack(blocks).tocsr()
inv = np.empty_like(order)
inv[order] = np.arange(len(order))
M = M[inv]
assert M.shape == a.shape

o = a.obs.copy()
o["scrublet_doublet"] = f.loc[o.index, "scrublet_doublet"].values
o["scDblFinder_score"] = res.loc[o.index, "scDblFinder_score"].values
o["scDblFinder_class"] = res.loc[o.index, "scDblFinder_class"].values
o["consensus_doublet"] = o.scrublet_doublet & (o.scDblFinder_class == "doublet")
o["decontX_contamination"] = res.loc[o.index, "decontX_contamination"].values
o["umis_raw"] = np.asarray(a.X.sum(1)).ravel()
o["umis_decontx"] = np.asarray(M.sum(1)).ravel()
for c in E.columns:
    if c not in o.columns:
        o[c] = E[c].values
o["annotation_draft"] = o["annotation_draft"].astype(str)
o["dx_lt200"] = o.umis_decontx < 200

# per-sample summary
g = o.groupby("sample", observed=True)
ps = pd.DataFrame({"suspension": g.suspension.first(), "digestion": g.digestion.first(), "n_kept": g.size(),
                   "scrublet": g.scrublet_doublet.mean(), "scDblFinder": g.apply(lambda x: (x.scDblFinder_class == "doublet").mean()),
                   "consensus": g.consensus_doublet.mean(), "median_contamination": g.decontX_contamination.median(),
                   "umis_kept": g.umis_decontx.sum() / g.umis_raw.sum(), "dx_lt200": g.dx_lt200.sum()})

# Alb diagnostic
sym = a.var["symbol"].astype(str).values
alb = np.flatnonzero(sym == "Alb")[0]
cp = lambda X, u: np.asarray(X[:, alb].todense()).ravel() / np.maximum(u, 1) * 1e4  # noqa: E731
o["Alb_cp10k_raw"] = cp(a.X, o.umis_raw.values)
o["Alb_cp10k_decontx"] = cp(M, o.umis_decontx.values)
epi = o.annotation_draft.isin(["Hepatocyte", "Cholangiocyte"])
hep = o.annotation_draft == "Hepatocyte"
ps["nonepi_Alb_cp10k_raw"] = o[~epi].groupby("sample", observed=True).Alb_cp10k_raw.mean()
ps["nonepi_Alb_cp10k_dx"] = o[~epi].groupby("sample", observed=True).Alb_cp10k_decontx.mean()
alb_raw = np.asarray(a.X[:, alb].todense()).ravel()
alb_dx = np.asarray(M[:, alb].todense()).ravel()
ps["hep_Alb_kept"] = pd.Series(alb_dx[hep.values], o.index[hep.values]).groupby(o.loc[hep, "sample"].values).sum() / \
    pd.Series(alb_raw[hep.values], o.index[hep.values]).groupby(o.loc[hep, "sample"].values).sum()
ps["nonepi_Alb_removed"] = 1 - pd.Series(alb_dx[~epi.values]).groupby(o.loc[~epi, "sample"].values).sum() / \
    pd.Series(alb_raw[~epi.values]).groupby(o.loc[~epi, "sample"].values).sum()
ps.sort_values(["suspension", "digestion"]).to_csv(f"{D}/per_sample_summary.tsv", sep="\t", float_format="%.4g")
print(ps.sort_values(["suspension", "digestion"]).to_string(float_format=lambda v: f"{v:.3f}"), flush=True)

# per-cluster (r3) summary for cluster-level doublet / contamination review
g = o.groupby("leiden_totalvi_r3", observed=True)
pc = pd.DataFrame({"annotation_draft": g.annotation_draft.agg(lambda x: x.mode().iat[0]), "n": g.size(),
                   "frac_nucleus": g.apply(lambda x: (x.suspension == "nucleus").mean()),
                   "scDblFinder": g.apply(lambda x: (x.scDblFinder_class == "doublet").mean()),
                   "scDblFinder_score_median": g.scDblFinder_score.median(), "consensus": g.consensus_doublet.mean(),
                   "median_contamination": g.decontX_contamination.median(), "dx_lt200": g.dx_lt200.mean()})
pc = pc.sort_values("scDblFinder", ascending=False)
pc.to_csv(f"{D}/per_cluster_r3_summary.tsv", sep="\t", float_format="%.4g")
print(pc.head(25).to_string(float_format=lambda v: f"{v:.3f}"), flush=True)

tot = {"cells": len(o), "scrublet": int(o.scrublet_doublet.sum()), "scDblFinder": int((o.scDblFinder_class == "doublet").sum()),
       "consensus": int(o.consensus_doublet.sum()), "dx_lt200": int(o.dx_lt200.sum()),
       "umis_kept": float(o.umis_decontx.sum() / o.umis_raw.sum()),
       "nonepi_Alb_cp10k_raw": float(o.Alb_cp10k_raw[~epi].mean()), "nonepi_Alb_cp10k_dx": float(o.Alb_cp10k_decontx[~epi].mean())}
print(tot, flush=True)
o.to_parquet(f"{D}/cells_doublets_ambient.parquet")

# figures: contamination and doublet score by draft type, split by suspension
types = o.groupby("annotation_draft").size().sort_values(ascending=False).index
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
for ax, col, ttl in [(axs[0], "decontX_contamination", "decontX contamination"),
                     (axs[1], "scDblFinder_score", "scDblFinder score")]:
    for j, (susp, colr) in enumerate([("cell", "tab:blue"), ("nucleus", "tab:orange")]):
        data = [o.loc[(o.annotation_draft == t) & (o.suspension == susp), col].dropna().values for t in types]
        pos = np.arange(len(types)) + (j - 0.5) * 0.4
        keep = [i for i, d in enumerate(data) if len(d) > 10]
        vp = ax.violinplot([data[i] for i in keep], positions=pos[keep], widths=0.38, showmedians=True)
        for b in vp["bodies"]:
            b.set_color(colr)
    ax.set_xticks(range(len(types)), types, rotation=90)
    ax.set_title(f"{ttl} (blue cells, orange nuclei)")
x = np.arange(len(types))
axs[2].bar(x - 0.2, [o.Alb_cp10k_raw[o.annotation_draft == t].mean() for t in types], 0.4, label="raw")
axs[2].bar(x + 0.2, [o.Alb_cp10k_decontx[o.annotation_draft == t].mean() for t in types], 0.4, label="decontX")
axs[2].set_yscale("log")
axs[2].set_xticks(x, types, rotation=90)
axs[2].set_title("Alb CP10k by draft type")
axs[2].legend()
fig.tight_layout()
fig.savefig(f"{D}/figures/doublets_ambient_by_type.png", dpi=200)
fig.savefig(f"{D}/figures/doublets_ambient_by_type.pdf")

Em = ad.read_h5ad(f"{OUT}/03_integration/embeddings.h5ad")
U = Em.obsm["X_umap_totalvi"]
fig, axs = plt.subplots(1, 4, figsize=(24, 6))
for ax, (v, ttl, kw) in zip(axs, [(o.scDblFinder_score.values, "scDblFinder score", dict(cmap="viridis")),
                                  (o.consensus_doublet.values.astype(float), "consensus doublet", dict(cmap="Reds")),
                                  (o.decontX_contamination.values, "decontX contamination", dict(cmap="magma", vmax=0.6)),
                                  (np.log1p(o.Alb_cp10k_decontx.values), "Alb log1p CP10k after decontX", dict(cmap="viridis"))]):
    ordr = np.argsort(v)
    sc_ = ax.scatter(U[ordr, 0], U[ordr, 1], c=v[ordr], s=0.2, rasterized=True, **kw)
    ax.set_title(ttl)
    ax.axis("off")
    fig.colorbar(sc_, ax=ax, shrink=0.6)
fig.savefig(f"{D}/figures/umap_doublets_ambient.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{D}/figures/umap_doublets_ambient.pdf", dpi=300, bbox_inches="tight")

h = ad.AnnData(X=M, obs=o, var=a.var.copy(), layers={"counts": a.X}, obsm={k: v for k, v in a.obsm.items()})
h.obsm["X_totalvi_draft"] = Em.obsm["X_totalvi"]
h.obsm["X_umap_totalvi_draft"] = U
h.uns["note"] = ("DRAFT 2026-11-11. X = decontX counts (celda, per sample, z = TotalVI Leiden r1, defaults, varGenes 5000, seed 12345); "
                 "layers counts = raw UMI. Cells = draft QC (authors 3-MAD rule). consensus_doublet = Scrublet & scDblFinder (not removed). "
                 "var_names = Ensembl IDs (Ens93), var symbol. obsm adt = 173 ADT raw counts (see adt_measured).")
h.write_h5ad(f"{D}/guilliams2022_fig1b_decontx.h5ad.tmp", compression="gzip")
os.replace(f"{D}/guilliams2022_fig1b_decontx.h5ad.tmp", f"{D}/guilliams2022_fig1b_decontx.h5ad")
print("done", h.shape, flush=True)
