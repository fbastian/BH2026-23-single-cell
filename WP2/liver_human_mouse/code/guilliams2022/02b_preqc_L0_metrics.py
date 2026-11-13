"""Step 02b - QC v1 inputs (MASLD final-QC recipe ported to mouse): pre-QC lineage L0 per cell + QC metrics.
Nothing is removed here.

L0 (as the MASLD per-study L0, marker-only): per sample, log-normalise, 2000 HVG, 30 PCs, Leiden r=1; each cluster
gets the lineage with the highest mean marker score (scanpy score_genes on the sample). Margin = best - second best
score, kept for review. Lineages: immune, epithelial (hepatocyte + cholangiocyte), endothelial, stromal.

Metrics: n_counts, n_genes, log1p of both, pct_mito (13 protein-coding mt- genes), pct_top20, pct_hb (haemoglobin
chains), pct_ribo (Rps/Rpl), complexity = log10 genes / log10 counts.

Input : <OUT>/01_raw/guilliams2022_fig1b_raw.h5ad
Output: <OUT>/02_qc_v1/cells_preqc.parquet (all 389,056 cells: obs + metrics + L0 + cluster + scores)
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/02_qc_v1"
os.makedirs(D, exist_ok=True)

MARKERS = {
    "immune": ["Ptprc", "Cd52", "Coro1a", "Laptm5", "Lcp1", "Cd53"],
    "epithelial": ["Alb", "Apoa2", "Ttr", "Serpina1a", "Epcam", "Krt19", "Krt7"],
    "endothelial": ["Kdr", "Stab2", "Ptprb", "Cdh5", "Clec4g", "Pecam1"],
    "stromal": ["Dcn", "Reln", "Col1a1", "Col3a1", "Lrat", "Pdgfrb", "Msln", "Upk3b"],
}
MT13 = ["mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4",
        "mt-Nd5", "mt-Nd6", "mt-Cytb"]
HB = ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Hbb-b1", "Hbb-b2"]

a = ad.read_h5ad(f"{OUT}/01_raw/guilliams2022_fig1b_raw.h5ad")
sym = a.var["symbol"].values
a.var_names = pd.Index(a.var["symbol"].astype(str)).str.cat(a.var_names, sep="|")  # unique keys
s2i = {s: i for i, s in enumerate(sym)}  # first occurrence
missing = [g for gs in MARKERS.values() for g in gs if g not in s2i] + [g for g in MT13 if g not in s2i]
print("missing marker / mt genes:", missing, flush=True)


def frac(X, genes, tot):
    idx = [s2i[g] for g in genes if g in s2i]
    return 100 * np.asarray(X[:, idx].sum(1)).ravel() / np.maximum(tot, 1)


X = a.X
tot = np.asarray(X.sum(1)).ravel()
m = a.obs[["sample", "gsm", "digestion", "suspension", "chemistry", "cite_seq", "condition", "mouse"]].copy()
m["n_counts"] = tot
m["n_genes"] = np.diff(X.indptr)
m["log1p_n_counts"] = np.log1p(m["n_counts"])
m["log1p_n_genes"] = np.log1p(m["n_genes"])
m["pct_mito"] = frac(X, MT13, tot)
m["pct_hb"] = frac(X, HB, tot)
m["pct_ribo"] = frac(X, [s for s in s2i if s.startswith(("Rps", "Rpl")) and "-" not in s], tot)
m["complexity"] = np.log10(np.maximum(m["n_genes"], 1)) / np.log10(np.maximum(tot, 2))
# top-20 share per cell
top20 = np.zeros(a.n_obs)
for i in range(a.n_obs):
    d = X.data[X.indptr[i]:X.indptr[i + 1]]
    top20[i] = d[np.argpartition(d, -20)[-20:]].sum() if len(d) > 20 else d.sum()
m["pct_top20"] = 100 * top20 / np.maximum(tot, 1)
print("metrics done", flush=True)

lab = pd.Series(index=a.obs_names, dtype=object)
clu = pd.Series(index=a.obs_names, dtype=object)
marg = pd.Series(index=a.obs_names, dtype=float)
scores = pd.DataFrame(index=a.obs_names, columns=list(MARKERS), dtype=float)
for s, idx in a.obs.groupby("sample", observed=True).indices.items():
    b = a[idx].copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    for L, gs in MARKERS.items():
        keys = [k for k in b.var_names if k.split("|")[0] in gs]
        sc.tl.score_genes(b, keys, score_name=f"sc_{L}", random_state=0)
    if b.n_obs >= 50:
        sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat")
        c = b[:, b.var["highly_variable"]].copy()
        sc.pp.scale(c, max_value=10)
        sc.pp.pca(c, n_comps=min(30, c.n_obs - 1), random_state=0)
        sc.pp.neighbors(c, random_state=0)
        sc.tl.leiden(c, resolution=1.0, random_state=0, flavor="igraph", n_iterations=2)
        cl = c.obs["leiden"].astype(str).values
    else:
        cl = np.array(["0"] * b.n_obs)
    S = b.obs[[f"sc_{L}" for L in MARKERS]].copy()
    S.columns = list(MARKERS)
    S["cl"] = cl
    means = S.groupby("cl").mean()
    best = means.idxmax(axis=1)
    srt = np.sort(means.values, axis=1)
    margin = pd.Series(srt[:, -1] - srt[:, -2], index=means.index)
    lab.iloc[idx] = best.loc[cl].values
    clu.iloc[idx] = [f"{s}|{c_}" for c_ in cl]
    marg.iloc[idx] = margin.loc[cl].values
    scores.iloc[idx] = S[list(MARKERS)].values
    print(f"{s:8s} n={b.n_obs:6d} clusters={len(means):3d} " +
          " ".join(f"{k}={v}" for k, v in pd.Series(best.loc[cl].values).value_counts().items()), flush=True)

m["L0"] = lab.astype("category")
m["L0_cluster"] = clu.values
m["L0_margin"] = marg.values
for L in MARKERS:
    m[f"score_{L}"] = scores[L].astype(float).values
m.to_parquet(f"{D}/cells_preqc.parquet")
print(pd.crosstab(m["digestion"], m["L0"]).to_string(), flush=True)
print("written", m.shape, flush=True)
