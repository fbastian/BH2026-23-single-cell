"""Step 12 - check the cDC2 <-> pDC streak on the strict TotalVI UMAP (user, 2026-11-12).
Zoom on the DC region: r3 clusters, labels, doublet score, sample/digestion/gate, markers for pDC, cDC2, cDC1,
transitional DC (tDC: Siglech+ Cx3cr1+ Zeb2+ Itgax+, Leylek 2019 / Sulczewski 2023), proliferation and doublet partners.
Output <OUT>/11_totalvi_strict/dc_streak/: umap_dc_zoom.png, r3_dc_clusters.tsv."""
import os
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/11_totalvi_strict/dc_streak"; os.makedirs(D, exist_ok=True)
E = ad.read_h5ad(f"{OUT}/11_totalvi_strict/embeddings.h5ad")
U = E.obsm["X_umap_totalvi"]; lab = E.obs.annotation_curated.astype(str)
dc = lab.isin(["pDCs", "cDC2", "cDC1", "Mig. cDCs"]).values
lo, hi = np.percentile(U[dc], 0.5, 0), np.percentile(U[dc], 99.5, 0)
pad = (hi - lo) * 0.08; lo, hi = lo - pad, hi + pad
box = ((U >= lo) & (U <= hi)).all(1)
print("cells in DC box", box.sum(), flush=True)
A = ad.read_h5ad(f"{OUT}/10_strict/guilliams2022_fig1b_strict.h5ad", backed="r")
assert (A.obs_names == E.obs_names).all()
m = A[np.flatnonzero(box)].to_memory(); s = ad.AnnData(X=m.X.copy(), obs=m.obs[[]].copy(), var=A.var.copy()); del m
s.obs = E.obs.loc[s.obs_names].copy(); s.obsm["X_umap"] = U[box]
sc.pp.normalize_total(s, target_sum=1e4); sc.pp.log1p(s)
s.var_names = s.var["symbol"].astype(str).values; s.var_names_make_unique()
G = {"pDC": ["Siglech", "Ccr9", "Bst2", "Ly6d", "Klk1", "Cox6a2", "Tcf4"],
     "cDC2": ["Cd209a", "Sirpa", "Itgam", "Clec4a4", "Mgl2"], "cDC1": ["Xcr1", "Clec9a"],
     "tDC/preDC": ["Cx3cr1", "Zeb2", "Ly6c2", "Itgax", "Lyz2", "Mki67"],
     "partner": ["Cd3e", "Cd79a", "Clec4f", "Alb", "Stab2", "Ptprc"]}
genes = [g for v in G.values() for g in v if g in s.var_names]
s.obs["log10_umis"] = np.log10(s.obs.umis_final.astype(float) + 1)
cols = ["annotation_curated", "leiden_clean_r3", "scDblFinder_score", "log10_umis", "digestion", "sort_gate", "sample", "bridge_type"] + genes
fig = sc.pl.umap(s, color=cols, ncols=6, s=4, show=False, return_fig=True, legend_fontsize=5, frameon=False, legend_loc="on data")
fig.savefig(f"{D}/umap_dc_zoom.png", dpi=150, bbox_inches="tight")
r3 = s.obs.leiden_clean_r3.astype(str)
X = pd.DataFrame(s[:, genes].X.toarray(), index=s.obs_names, columns=genes)
t = X.groupby(r3.values).mean()
t.insert(0, "n", r3.value_counts().reindex(t.index).values)
t.insert(1, "label", s.obs.groupby(r3.values).annotation_curated.agg(lambda x: x.astype(str).mode().iat[0]))
t.insert(2, "frac_label", s.obs.groupby(r3.values).annotation_curated.agg(lambda x: (x.astype(str) == x.astype(str).mode().iat[0]).mean()))
t.insert(3, "dbl_score", s.obs.groupby(r3.values).scDblFinder_score.median())
t.insert(4, "umis", s.obs.groupby(r3.values).umis_final.median())
t.insert(5, "exvivo", s.obs.groupby(r3.values).digestion.agg(lambda x: (x.astype(str) == "Ex Vivo").mean()))
t.insert(6, "top_sample", s.obs.groupby(r3.values)["sample"].agg(lambda x: f"{x.astype(str).mode().iat[0]} {(x.astype(str)==x.astype(str).mode().iat[0]).mean():.2f}"))
t.insert(7, "umap", [f"{a:.1f},{b:.1f}" for a, b in s.obs.groupby(r3.values).apply(lambda x: s.obsm["X_umap"][s.obs_names.get_indexer(x.index)].mean(0))])
t = t[t.n >= 100].sort_values("label")
t.to_csv(f"{D}/r3_dc_clusters.tsv", sep="\t", float_format="%.3g")
pd.set_option("display.width", 400)
print(t.to_string(float_format=lambda v: f"{v:.2f}"), flush=True)
