"""Step 05 - DRAFT annotation of the TotalVI clusters with the paper's Table S1 (mmc1) mouse DEGs (17 Fig 1B types).

Cluster-level scoring (fast, for the poster):
  1. Leiden r1/r2/r3 on X_totalvi (GPU).
  2. Mean log-normalised expression per r3 cluster (pseudobulk).
  3. Per type: top 50 Table S1 DEGs by the authors' 'score' column; per gene z-score across clusters; type score =
     mean z over the markers found.
  4. Label = best type; margin = best - second best; flag 'low_margin' when margin < 0.5.
Outputs <OUT>/05_annotation_draft/: cluster_scores.tsv, cluster_labels.tsv, cells_annotation_draft.parquet,
figures/umap_annotation_draft.{png,pdf}, figures/cluster_type_heatmap.{png,pdf}.
"""
import os
import sys

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
S1 = "/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/1-s2.0-S0092867421014811-mmc1.xlsx"
D = f"{OUT}/05_annotation_draft"
os.makedirs(f"{D}/figures", exist_ok=True)
NTOP = 50

E = ad.read_h5ad(f"{OUT}/03_integration/embeddings.h5ad")
import rapids_singlecell as rsc

t = ad.AnnData(obs=E.obs[[]].copy(), obsm={"X_totalvi": E.obsm["X_totalvi"].astype(np.float32)})
rsc.pp.neighbors(t, use_rep="X_totalvi")
for r in (2.0, 3.0):
    rsc.tl.leiden(t, resolution=r, random_state=0, key_added=f"r{r:g}")
    E.obs[f"leiden_totalvi_r{r:g}"] = t.obs[f"r{r:g}"].astype(str).values
print("clusters r1/r2/r3:", [E.obs[f"leiden_totalvi_r{r}"].nunique() for r in (1, 2, 3)], flush=True)

a = ad.read_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad")
assert (a.obs_names == E.obs_names).all()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sym = a.var["symbol"].astype(str).values

# markers: top NTOP DEGs per mouse type, mapped to columns (all Ensembl IDs sharing the symbol)
x = pd.read_excel(S1, sheet_name=None)
markers = {}
for k, v in x.items():
    if k.startswith("Mouse") and k.endswith("DEGs"):
        typ = k.replace("Mouse ", "").replace(" DEGs", "")
        g = v.sort_values("score", ascending=False).iloc[:NTOP, 0].astype(str)
        markers[typ] = [s for s in g if s in set(sym)]
print({k: len(v) for k, v in markers.items()}, flush=True)

key = "leiden_totalvi_r3"
cl = E.obs[key].astype(str).values
cats = sorted(set(cl), key=int)
ind = sp.csr_matrix((np.ones(len(cl)), (pd.Categorical(cl, categories=cats).codes, np.arange(len(cl)))),
                    shape=(len(cats), len(cl)))
n = np.asarray(ind.sum(1)).ravel()
need = sorted({s for v in markers.values() for s in v})
cols = np.where(np.isin(sym, need))[0]
M = np.asarray((ind @ a.X[:, cols]).todense()) / n[:, None]  # cluster x gene mean log-norm
G = pd.DataFrame(M, index=cats, columns=sym[cols]).T.groupby(level=0).mean().T  # collapse duplicate symbols
Z = (G - G.mean()) / G.std().replace(0, np.nan)
score = pd.DataFrame({t_: Z[g].mean(1) for t_, g in markers.items()})
score.to_csv(f"{D}/cluster_scores.tsv", sep="\t")
srt = np.sort(score.values, 1)
lab = pd.DataFrame({"n_cells": n.astype(int), "label": score.idxmax(1), "score": srt[:, -1],
                    "second": score.apply(lambda r: r.drop(r.idxmax()).idxmax(), axis=1),
                    "margin": srt[:, -1] - srt[:, -2]}, index=cats)
lab["low_margin"] = lab["margin"] < 0.5
# composition per cluster for context
for c in ("suspension", "digestion"):
    tab = pd.crosstab(E.obs[key], E.obs[c], normalize="index")
    for k in tab.columns:
        lab[f"frac_{k}"] = tab.loc[lab.index, k].round(3).values
lab.to_csv(f"{D}/cluster_labels.tsv", sep="\t")
print(lab.to_string(), flush=True)

E.obs["annotation_draft"] = pd.Categorical(lab.loc[cl, "label"].values)
E.obs["annotation_draft_low_margin"] = lab.loc[cl, "low_margin"].values
E.obs[["leiden_totalvi_r1", "leiden_totalvi_r2", "leiden_totalvi_r3", "annotation_draft",
       "annotation_draft_low_margin", "suspension", "digestion", "sample"]].to_parquet(
    f"{D}/cells_annotation_draft.parquet")
comp = E.obs.groupby(["annotation_draft", "suspension"], observed=True).size().unstack(fill_value=0)
comp["total"] = comp.sum(1)
comp.sort_values("total", ascending=False).to_csv(f"{D}/composition_draft.tsv", sep="\t")
print(comp.sort_values("total", ascending=False).to_string(), flush=True)

# figures
p = ad.AnnData(obs=E.obs.copy(), obsm={"X_umap": E.obsm["X_umap_totalvi"]})
fig, axs = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(p, color="annotation_draft", ax=axs[0], show=False, legend_loc="on data", legend_fontsize=6,
           frameon=False, title="TotalVI: draft annotation (Table S1 markers)", s=1)
sc.pl.umap(p, color="annotation_draft", ax=axs[1], show=False, frameon=False, title="", s=1)
fig.savefig(f"{D}/figures/umap_annotation_draft.png", dpi=200, bbox_inches="tight")
for ax in axs:
    for c in ax.collections:
        c.set_rasterized(True)
fig.savefig(f"{D}/figures/umap_annotation_draft.pdf", dpi=300, bbox_inches="tight")
order = lab.sort_values(["label", "score"], ascending=[True, False]).index
fig, ax = plt.subplots(figsize=(8, max(6, len(order) * 0.12)))
im = ax.imshow(score.loc[order].values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_yticks(range(len(order)), [f"{c} {lab.loc[c, 'label']}" for c in order], fontsize=4)
ax.set_xticks(range(score.shape[1]), score.columns, rotation=90, fontsize=6)
fig.colorbar(im, ax=ax, label="mean marker z")
fig.savefig(f"{D}/figures/cluster_type_heatmap.png", dpi=200, bbox_inches="tight")
fig.savefig(f"{D}/figures/cluster_type_heatmap.pdf", bbox_inches="tight")
print("done", flush=True)
