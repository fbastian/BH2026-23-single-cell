"""Step 09 - why streaks? Per cell: bridge score = fraction of 15 nearest neighbours (TotalVI clean latent) in a different
Leiden r1 cluster. Cross-tab against scDblFinder, low UMI, decontX lt200, CS140, suspension. Then re-draw UMAP after removing
(a) per-cell scDblFinder doublets, (b) + bridge cells (score > 0.5) - an approximation of the authors' strict QC.
Output 08_totalvi_clean/streaks/: tables + figures."""
import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cuml.neighbors import NearestNeighbors

D = "/mnt/masld_vault/cross_species/guilliams2022/08_totalvi_clean"
O = f"{D}/streaks"
import os
os.makedirs(O, exist_ok=True)
E = ad.read_h5ad(f"{D}/embeddings.h5ad")
Z = E.obsm["X_totalvi"].astype(np.float32)
nn = NearestNeighbors(n_neighbors=16).fit(Z)
_, idx = nn.kneighbors(Z)
idx = np.asarray(idx)[:, 1:]
cl = E.obs["leiden_clean_r1"].astype(str).values
E.obs["bridge"] = (cl[idx] != cl[:, None]).mean(1)
o = E.obs
o["low_umi"] = o["umis_final"] < 1000
o["bridge_cell"] = o["bridge"] > 0.5
flags = ["scDblFinder_doublet", "dx_lt200", "low_umi", "ambient_dominated"]
rows = []
for f in flags + ["bridge_cell"]:
    v = o[f].astype(bool)
    rows.append(dict(flag=f, n=int(v.sum()), pct_all=v.mean() * 100,
                     pct_among_bridge=v[o.bridge_cell].mean() * 100, pct_among_nonbridge=v[~o.bridge_cell].mean() * 100))
t = pd.DataFrame(rows)
t.to_csv(f"{O}/flags_vs_bridge.tsv", sep="\t", index=False, float_format="%.2f")
print(t.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
print("bridge score by scDblFinder class:\n", o.groupby("scDblFinder_doublet").bridge.describe().round(3))
print("bridge by suspension:\n", o.groupby("suspension", observed=True).bridge_cell.mean().round(3))

U = E.obsm["X_umap_totalvi"]
fig, axs = plt.subplots(1, 4, figsize=(28, 7))
for ax, (v, ttl, cm) in zip(axs, [(o.bridge.values, "bridge score (frac. neighbours in other cluster)", "magma"),
                                  (o.scDblFinder_doublet.values.astype(float), "scDblFinder doublet (per cell)", "Reds"),
                                  (np.log10(o.umis_final.values + 1), "log10 UMI (cleaned)", "viridis"),
                                  (o.decontX_contamination_final.fillna(0).values, "decontX contamination (tier-1 samples)", "magma")]):
    r = np.argsort(v)
    s = ax.scatter(U[r, 0], U[r, 1], c=v[r], s=0.2, cmap=cm, rasterized=True)
    ax.set_title(ttl); ax.axis("off"); fig.colorbar(s, ax=ax, shrink=0.6)
fig.savefig(f"{O}/umap_streak_drivers.png", dpi=150, bbox_inches="tight")

def redraw(mask, tag):
    g = ad.AnnData(obs=o.loc[mask, ["annotation_curated"]].copy(), obsm={"Z": Z[mask]})
    rsc.pp.neighbors(g, use_rep="Z")
    rsc.tl.umap(g, random_state=0)
    return g

ga = redraw(~o.scDblFinder_doublet.values, "a")
gb = redraw(~(o.scDblFinder_doublet.values | o.bridge_cell.values | o.dx_lt200.values), "b")
fig, axs = plt.subplots(1, 3, figsize=(30, 9))
full = ad.AnnData(obs=o[["annotation_curated"]].copy(), obsm={"X_umap": U})
for ax, (g, ttl) in zip(axs, [(full, f"all clean cells ({len(o):,})"), (ga, f"- scDblFinder doublets ({ga.n_obs:,})"),
                              (gb, f"- doublets, bridge cells, <200 UMI after decontX ({gb.n_obs:,})")]):
    sc.pl.umap(g, color="annotation_curated", ax=ax, show=False, frameon=False, s=1, legend_loc="on data",
               legend_fontsize=6, title=ttl)
    for c in ax.collections:
        c.set_rasterized(True)
fig.savefig(f"{O}/umap_strict_qc_approx.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{O}/umap_strict_qc_approx.pdf", dpi=300, bbox_inches="tight")
print("cells after a:", ga.n_obs, "after b:", gb.n_obs)
comp = pd.concat([o.annotation_curated.value_counts().rename("clean"),
                  gb.obs.annotation_curated.value_counts().rename("strict_approx")], axis=1)
comp.to_csv(f"{O}/composition_strict_approx.tsv", sep="\t")
print(comp.to_string())
