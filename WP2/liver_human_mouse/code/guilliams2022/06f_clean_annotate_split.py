"""Step 06f - clean object, re-annotation, split into single cells and single nuclei (user go 2026-11-11).

Ambient (from the 06d ladder):
  - decontX (1, 300) for the 13 in vivo cell samples with non-epithelial Alb > 5 CP10k (06e outputs);
  - raw counts for ex vivo, the clean in vivo samples and all nuclei (decontX removed hepatocyte signal in nuclei);
  - CS140: raw counts + flag ambient_dominated (ambient Alb 440 CP10k; decontX cost hepatocytes ~30 % Alb at every strength;
    decision open, default = MASLD raw_by_decision precedent).
  - dx_lt200: < 200 UMI after decontX (decontX samples only, neutrophils exempt) - flagged, not removed.
Doublets: the draft QC already removed Scrublet calls, so the MASLD consensus rule cannot apply. Removed = the 8 TotalVI
Leiden r3 clusters with > 65 % scDblFinder doublets (authors' approach: doublet clusters removed by hand). Per-cell scDblFinder
calls elsewhere are flagged only.
Annotation: Table S1 (mmc1) top-50 mouse DEGs, cluster-level (Leiden r3), log-normalised cleaned counts; z-scores computed
separately within cell-majority and nucleus-majority clusters (mixing them mislabelled hepatocyte nuclei in the draft).
Outputs <OUT>/07_clean/: guilliams2022_fig1b_clean.h5ad (all), _cells.h5ad, _nuclei.h5ad, cluster_labels_clean.tsv,
composition_clean.tsv, label_changes.tsv, figures/.
"""
import os

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/mnt/masld_vault/cross_species/guilliams2022"
DX = f"{OUT}/06_doublets_ambient/final_dx"
D = f"{OUT}/07_clean"
os.makedirs(f"{D}/figures", exist_ok=True)
S1 = "/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/1-s2.0-S0092867421014811-mmc1.xlsx"
TIER1 = "CS114 CS137 CS138 CS139 CS144 CS48 CS52 CS88 CS89 CS92 CS93 CS96 CS97".split()
DOUBLET_CLUSTERS = ["45", "41", "50", "12", "5", "71", "23", "39"]
NTOP = 50

a = ad.read_h5ad(f"{OUT}/06_doublets_ambient/guilliams2022_fig1b_decontx.h5ad")
raw = a.layers["counts"].tocsr()
obs = a.obs.copy()
obs["decontx_setting"] = "raw"
obs["decontX_contamination_final"] = np.nan
rows, blocks = [], []
for s in TIER1:
    idx = np.flatnonzero(obs["sample"].astype(str).values == s)
    m = sp.load_npz(f"{DX}/{s}.dx300.npz").tocsr()
    c = pd.read_parquet(f"{DX}/{s}.dx300.parquet")
    assert m.shape[0] == len(idx) and (c.index == obs.index[idx]).all(), s
    rows.append(idx)
    blocks.append(m)
    obs.iloc[idx, obs.columns.get_loc("decontx_setting")] = "decontX_1_300"
    obs.iloc[idx, obs.columns.get_loc("decontX_contamination_final")] = c.iloc[:, 0].values
# assemble final counts: raw rows except tier-1 rows replaced
repl = np.concatenate(rows)
keep_raw = np.setdiff1d(np.arange(a.n_obs), repl)
parts = sp.vstack([raw[keep_raw].astype(np.float32), sp.vstack(blocks).astype(np.float32)]).tocsr()
order = np.concatenate([keep_raw, repl])
inv = np.empty_like(order)
inv[order] = np.arange(len(order))
X = parts[inv]
obs["umis_final"] = np.asarray(X.sum(1)).ravel()
obs["ambient_dominated"] = obs["sample"].astype(str) == "CS140"
obs["dx_lt200"] = (obs.decontx_setting != "raw") & (obs.umis_final < 200) & (obs.annotation_draft != "Neutrophil")
obs["scDblFinder_doublet"] = obs.scDblFinder_class == "doublet"
obs["doublet_cluster"] = obs.leiden_totalvi_r3.astype(str).isin(DOUBLET_CLUSTERS)
print("cells", len(obs), "doublet-cluster removed", int(obs.doublet_cluster.sum()), "dx_lt200 flagged", int(obs.dx_lt200.sum()),
      "CS140 flagged", int(obs.ambient_dominated.sum()), flush=True)

c = ad.AnnData(X=X, obs=obs, var=a.var.copy(), layers={"counts": raw}, obsm=dict(a.obsm))
c = c[~c.obs.doublet_cluster.values].copy()
del a, raw, X, parts

# re-annotation on cleaned counts
n = c.copy()
sc.pp.normalize_total(n, target_sum=1e4)
sc.pp.log1p(n)
sym = n.var["symbol"].astype(str).values
x = pd.read_excel(S1, sheet_name=None)
markers = {}
for k, v in x.items():
    if k.startswith("Mouse") and k.endswith("DEGs"):
        g = v.sort_values("score", ascending=False).iloc[:NTOP, 0].astype(str)
        markers[k.replace("Mouse ", "").replace(" DEGs", "")] = [s for s in g if s in set(sym)]
cl = n.obs.leiden_totalvi_r3.astype(str).values
cats = sorted(set(cl), key=int)
codes = pd.Categorical(cl, categories=cats).codes
ind = sp.csr_matrix((np.ones(len(cl)), (codes, np.arange(len(cl)))), shape=(len(cats), len(cl)))
ncl = np.asarray(ind.sum(1)).ravel()
need = sorted({s for v in markers.values() for s in v})
cols = np.where(np.isin(sym, need))[0]
G = pd.DataFrame(np.asarray((ind @ n.X[:, cols]).todense()) / ncl[:, None], index=cats, columns=sym[cols]).T.groupby(level=0).mean().T
fnuc = pd.Series((n.obs.suspension.astype(str).values == "nucleus"), index=n.obs_names).groupby(cl).mean().reindex(cats)
grp = np.where(fnuc.values > 0.5, "nucleus", "cell")
score = pd.DataFrame(index=cats, columns=list(markers), dtype=float)
for gname in ("cell", "nucleus"):
    sel = grp == gname
    Gs = G.loc[sel]
    Z = (Gs - Gs.mean()) / Gs.std().replace(0, np.nan)
    for t_, gl in markers.items():
        score.loc[Gs.index, t_] = Z[gl].mean(1)
score.to_csv(f"{D}/cluster_scores_clean.tsv", sep="\t", float_format="%.4g")
srt = np.sort(score.values, 1)
lab = pd.DataFrame({"n_cells": ncl.astype(int), "z_group": grp, "frac_nucleus": fnuc.values, "label": score.idxmax(1),
                    "score": srt[:, -1], "second": score.apply(lambda r: r.drop(r.idxmax()).idxmax(), axis=1),
                    "margin": srt[:, -1] - srt[:, -2]}, index=cats)
lab["low_margin"] = lab.margin < 0.5
draft = c.obs.groupby(c.obs.leiden_totalvi_r3.astype(str), observed=True).annotation_draft.agg(lambda s: s.mode().iat[0])
lab["label_draft"] = draft.reindex(cats).values
lab.to_csv(f"{D}/cluster_labels_clean.tsv", sep="\t", float_format="%.4g")
lab[lab.label != lab.label_draft].to_csv(f"{D}/label_changes.tsv", sep="\t", float_format="%.4g")
print(lab.to_string(float_format=lambda v: f"{v:.2f}"), flush=True)
print("label changes draft -> clean:\n", lab[lab.label != lab.label_draft][["n_cells", "frac_nucleus", "label_draft", "label",
                                                                              "margin"]].to_string(), flush=True)
del n

c.obs["annotation_clean"] = pd.Categorical(lab.loc[cl, "label"].values)
c.obs["annotation_clean_low_margin"] = lab.loc[cl, "low_margin"].values
comp = c.obs.groupby(["annotation_clean", "suspension"], observed=True).size().unstack(fill_value=0)
comp["total"] = comp.sum(1)
comp = comp.sort_values("total", ascending=False)
comp.to_csv(f"{D}/composition_clean.tsv", sep="\t")
print(comp.to_string(), "\ntotal", int(comp.total.sum()), flush=True)
for k in list(c.obs.columns):
    if c.obs[k].dtype == object:
        c.obs[k] = c.obs[k].astype(str)
c.uns["note"] = ("DRAFT-CLEAN 2026-11-11. X = final counts (decontX (1,300) for 13 in vivo cell samples, raw otherwise; "
                 "obs decontx_setting); layers counts = raw UMI. Removed: 8 doublet clusters (TotalVI Leiden r3, >65 % scDblFinder). "
                 "Flags kept: scDblFinder_doublet, dx_lt200, ambient_dominated (CS140, raw). annotation_clean = Table S1 top-50 DEGs, "
                 "cluster level, z within cell / nucleus clusters. Embedding X_totalvi_draft = TotalVI on uncleaned counts.")
c.write_h5ad(f"{D}/guilliams2022_fig1b_clean.h5ad", compression="gzip")
print("written clean", c.shape, flush=True)

import rapids_singlecell as rsc  # noqa: E402

for susp, tag in (("cell", "cells"), ("nucleus", "nuclei")):
    s = c[c.obs.suspension.astype(str).values == susp].copy()
    t = ad.AnnData(obs=s.obs[[]].copy(), obsm={"X": s.obsm["X_totalvi_draft"].astype(np.float32)})
    rsc.pp.neighbors(t, use_rep="X")
    rsc.tl.umap(t, random_state=0)
    s.obsm["X_umap"] = np.asarray(t.obsm["X_umap"])
    if susp == "nucleus":
        for k in ("adt", "adt_measured"):
            s.obsm.pop(k, None)
    s.write_h5ad(f"{D}/guilliams2022_fig1b_clean_{tag}.h5ad", compression="gzip")
    fig, axs = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(s, color="annotation_clean", ax=axs[0], show=False, legend_loc="on data", legend_fontsize=6, frameon=False,
               title=f"{tag}: {s.n_obs:,} (draft TotalVI latent, cleaned labels)", s=2)
    sc.pl.umap(s, color="digestion", ax=axs[1], show=False, frameon=False, s=2)
    fig.savefig(f"{D}/figures/umap_clean_{tag}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("written", tag, s.shape, flush=True)
print("done", flush=True)
