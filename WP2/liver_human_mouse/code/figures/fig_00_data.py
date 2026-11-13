"""Figure data (2026-11-12): small tables so the figure scripts run in seconds.
1. Unintegrated baseline for Fig 2: same SAMap input cells (hs 50,400 + mm 53,370), genes = MGI 1:1 orthologs present in both
   (mouse symbol), counts concatenated, normalise_total 1e4 + log1p, 3,000 HVG (seurat, no batch key), PCA 50, kNN 15, UMAP seed 0 (rapids-singlecell, GPU).
2. Fig 1 dot plot statistics: mean log-normalised expression and fraction > 0 per annotation_aor type for the Table S1 top-3 genes.
3. Per-cell tables: mouse AoR UMAP + labels; SAMap joint UMAP + labels.
Outputs <D>/data/."""
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
from common import D
I = "/mnt/masld_vault/cross_species/samap_run/input"
# 1 unintegrated
mgi = pd.read_csv("/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/refs/orthologs/HOM_MouseHumanSequence.rpt", sep="\t", dtype=str)
one = mgi.groupby("DB Class Key").filter(lambda x: (x["NCBI Taxon ID"] == "10090").sum() == 1 and (x["NCBI Taxon ID"] == "9606").sum() == 1)
p = one.pivot_table(index="DB Class Key", columns="NCBI Taxon ID", values="Symbol", aggfunc="first").rename(columns={"9606": "hs", "10090": "mm"})
h = ad.read_h5ad(f"{I}/hs_samap_ready.h5ad"); m = ad.read_h5ad(f"{I}/mm_samap_ready.h5ad")
p = p[p.hs.isin(h.var_names) & p.mm.isin(m.var_names)]
print("1:1 orthologs in both objects:", len(p), flush=True)
hh = h[:, p.hs.values].copy(); hh.var_names = p.mm.values
mm_ = m[:, p.mm.values].copy()
hh.obs = pd.DataFrame({"species": "human", "label": h.obs.author_cell_type.astype(str).values}, index="hs_" + h.obs_names)
mm_.obs = pd.DataFrame({"species": "mouse", "label": m.obs.annotation_aor.astype(str).values}, index="mm_" + m.obs_names)
u = ad.concat([hh, mm_]); u.X = sp.csr_matrix(u.X)
sc.pp.normalize_total(u, target_sum=1e4); sc.pp.log1p(u)
sc.pp.highly_variable_genes(u, n_top_genes=3000, flavor="seurat")
sc.pp.pca(u, n_comps=50, mask_var="highly_variable", random_state=0)
import rapids_singlecell as rsc; rsc.pp.neighbors(u, n_neighbors=15, use_rep="X_pca"); rsc.tl.umap(u, random_state=0)
pd.DataFrame({"species": u.obs.species.values, "label": u.obs.label.values, "u1": u.obsm["X_umap"][:, 0], "u2": u.obsm["X_umap"][:, 1]},
             index=u.obs_names).to_parquet(f"{D}/data/unintegrated_umap.parquet")
print("unintegrated done", flush=True)
del u, hh, mm_, h, m
# 3a SAMap joint
j = ad.read_h5ad("/mnt/masld_vault/cross_species/samap_run/out/samap_joint.h5ad", backed="r")
sp_ = j.obs.species.astype(str).values
lab = np.where(sp_ == "hs", j.obs.hs_author_cell_type.astype(str).values, j.obs.mm_annotation_aor.astype(str).values)
pd.DataFrame({"species": np.where(sp_ == "hs", "human", "mouse"), "label": lab, "u1": j.obsm["X_umap"][:, 0], "u2": j.obsm["X_umap"][:, 1]},
             index=j.obs_names).to_parquet(f"{D}/data/samap_umap.parquet")
# 2 + 3b mouse
a = ad.read_h5ad("/mnt/masld_vault/cross_species/guilliams2022/13_annotation_of_record/guilliams2022_fig1b_aor.h5ad")
pd.DataFrame({"label": a.obs.annotation_aor.astype(str).values, "suspension": a.obs.suspension.astype(str).values,
              "digestion": a.obs.digestion.astype(str).values, "u1": a.obsm["X_umap"][:, 0], "u2": a.obsm["X_umap"][:, 1]},
             index=a.obs_names).to_parquet(f"{D}/data/mouse_aor_umap.parquet")
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.var_names = a.var.symbol.astype(str).values; a.var_names_make_unique()
S1 = "/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/1-s2.0-S0092867421014811-mmc1.xlsx"
x = pd.read_excel(S1, sheet_name=None); genes = {}
for k, v in x.items():
    if k.startswith("Mouse") and k.endswith("DEGs"):
        genes[k.replace("Mouse ", "").replace(" DEGs", "")] = [g for g in v.sort_values("score", ascending=False).iloc[:, 0].astype(str) if g in a.var_names][:3]
genes["Transitional DC"] = ["Ly6c2", "Cd209a", "Siglech"]; genes["Mesothelial"] = ["Msln", "Upk3b", "Lrrn4"]
pd.Series({k: ",".join(v) for k, v in genes.items()}).to_csv(f"{D}/data/dotplot_genes.tsv", sep="\t", header=False)
gl = list(dict.fromkeys(g for v in genes.values() for g in v))
lab = a.obs.annotation_aor.astype(str).values
X = a[:, gl].X; X = X.toarray() if sp.issparse(X) else X
mean = pd.DataFrame(X, columns=gl).groupby(lab).mean(); frac = pd.DataFrame(X > 0, columns=gl).groupby(lab).mean()
mean.to_csv(f"{D}/data/dotplot_mean.tsv", sep="\t"); frac.to_csv(f"{D}/data/dotplot_frac.tsv", sep="\t")
print("done", flush=True)
