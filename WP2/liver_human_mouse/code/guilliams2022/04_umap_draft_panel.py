"""Draft UMAP panel: rows unintegrated / Harmony (TotalVI row added when it exists); columns digestion, suspension,
CITE-seq, Leiden (Harmony) and marker expression. Writes <OUT>/03_integration/figures/umap_draft_panel.{pdf,png}"""
import os, sys
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/03_integration"; os.makedirs(f"{D}/figures", exist_ok=True)
E = ad.read_h5ad(f"{D}/embeddings_harmony.h5ad")
if os.path.exists(f"{D}/embeddings.h5ad"):
    T = ad.read_h5ad(f"{D}/embeddings.h5ad")
    if "X_umap_totalvi" in T.obsm: E.obsm["X_umap_totalvi"] = T.obsm["X_umap_totalvi"]
a = ad.read_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad")
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
MK = ["Alb", "Clec4f", "Ptprc", "Cd3e", "Cd79a", "Stab2", "Reln", "Krt19", "S100a8"]
sym = pd.Series(a.var_names, index=a.var["symbol"].values)
for g in MK:
    E.obs[g] = np.asarray(a[:, sym[g] if np.ndim(sym[g]) == 0 else sym[g].iloc[0]].X.todense()).ravel()
cols = ["digestion", "suspension", "cite_seq", "leiden_harmony_r1"] + MK
E.obs["cite_seq"] = E.obs["cite_seq"].astype(str).astype("category")
rows = [("X_umap_unintegrated", "Unintegrated"), ("X_umap_harmony", "Harmony")] + \
       ([("X_umap_totalvi", "TotalVI")] if "X_umap_totalvi" in E.obsm else [])
plt.rcParams.update({"font.size": 6, "pdf.fonttype": 42})
fig, axes = plt.subplots(len(rows), len(cols), figsize=(2.0 * len(cols), 2.1 * len(rows)))
idx = np.random.default_rng(0).permutation(E.n_obs)
for i, (key, name) in enumerate(rows):
    for j, c in enumerate(cols):
        ax = axes[i, j]
        sc.pl.embedding(E[idx], basis=key, color=c, ax=ax, show=False, size=0.5, frameon=False,
                        legend_loc="on data" if c.startswith("leiden") else "right margin",
                        legend_fontsize=4, title=f"{name}: {c}" if i == 0 or j == 0 else c,
                        color_map="viridis", vmax="p99" if c in MK else None)
fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(f"{D}/figures/umap_draft_panel.{ext}", dpi=200, bbox_inches="tight")
print("ok", E.n_obs)
