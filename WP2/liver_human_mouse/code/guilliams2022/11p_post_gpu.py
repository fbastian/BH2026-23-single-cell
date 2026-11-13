"""Step 11p - GPU neighbours/UMAP/Leiden r1-r3 on the clean TotalVI latent + UMAP panel (rapids-sc env)."""
import anndata as ad
import numpy as np
import scanpy as sc
import rapids_singlecell as rsc

OUT = "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/11_totalvi_strict"
E = ad.read_h5ad(f"{D}/embeddings.h5ad")
a = ad.read_h5ad(f"{OUT}/10_strict/guilliams2022_fig1b_strict.h5ad")
assert (E.obs_names == a.obs_names).all()
a.layers.pop("counts")
def log(m):
    print(m, flush=True)

g = ad.AnnData(obs=E.obs[[]].copy(), obsm={"X_totalvi": E.obsm["X_totalvi"]})
rsc.pp.neighbors(g, use_rep="X_totalvi")
rsc.tl.umap(g, random_state=0)
E.obsm["X_umap_totalvi"] = np.asarray(g.obsm["X_umap"])
for r in (1.0, 2.0, 3.0):
    rsc.tl.leiden(g, resolution=r, random_state=0, key_added="l")
    E.obs[f"leiden_clean_r{r:g}"] = g.obs["l"].astype(str).values
for k in E.obs.columns:
    if E.obs[k].dtype == object:
        E.obs[k] = E.obs[k].astype(str)
E.write_h5ad(f"{D}/embeddings.h5ad")
log(f"embeddings written; clusters r1/r2/r3 {[E.obs[f'leiden_clean_r{r}'].nunique() for r in (1, 2, 3)]}")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

n = a.copy()
sc.pp.normalize_total(n, target_sum=1e4)
sc.pp.log1p(n)
n.obsm["X_umap"] = E.obsm["X_umap_totalvi"]
n.obs = E.obs
n.var_names = n.var["symbol"].astype(str).values
n.var_names_make_unique()
cols = ["annotation_curated", "annotation_flag", "suspension", "digestion", "decontx_setting", "leiden_clean_r1",
        "Alb", "Clec4f", "Ptprc", "Stab2", "Reln", "Krt19", "Cd3e", "Cd79a", "S100a8", "Siglech"]
fig, axs = plt.subplots(4, 4, figsize=(28, 26))
for ax, c in zip(axs.ravel(), cols):
    sc.pl.umap(n, color=c, ax=ax, show=False, frameon=False, s=1,
               legend_loc="on data" if c in ("annotation_curated", "leiden_clean_r1") else "right margin", legend_fontsize=6)
    for coll in ax.collections:
        coll.set_rasterized(True)
fig.suptitle("TotalVI on strict-QC cells (10_strict) - curated draft labels", fontsize=14)
fig.savefig(f"{D}/figures/umap_totalvi_clean_panel.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{D}/figures/umap_totalvi_clean_panel.pdf", dpi=300, bbox_inches="tight")
log("finished")
