"""Step 03h - unintegrated PCA + Harmony baseline on GPU (rapids-singlecell), same HVG rule as 03_integrate.py.
Writes <OUT>/03_integration/embeddings_harmony.h5ad (obs + X_pca, X_pca_harmony, X_umap_unintegrated, X_umap_harmony,
leiden_harmony_r1). Run with CUDA_VISIBLE_DEVICES set to the free T4.
"""
import os
import sys
import time

import anndata as ad
import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import rmm
import scanpy as sc

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/03_integration"
os.makedirs(D, exist_ok=True)
rmm.reinitialize(managed_memory=False, pool_allocator=False)  # plain allocator (rsc Harmony trap, MASLD test17)
cp.cuda.set_allocator(rmm.allocators.cupy.rmm_cupy_allocator)
t0 = time.time()


def log(m):
    print(f"[{(time.time() - t0) / 60:6.1f} min] {m}", flush=True)


a = ad.read_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad")
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=4000, flavor="seurat", batch_key="sample")
a = a[:, a.var["highly_variable"]].copy()
log(f"HVG {a.shape}")
# scale + PCA on CPU: the dense scaled matrix (328k x 4000 float32, 5 GB) plus copies does not fit a T4
sc.pp.scale(a, max_value=10)
sc.pp.pca(a, n_comps=50, random_state=0)
log("PCA (CPU) done")
a = ad.AnnData(obs=a.obs.copy(), obsm={"X_pca": a.obsm["X_pca"].astype(np.float32)})
E = ad.AnnData(obs=a.obs.copy())
E.obsm["X_pca"] = np.asarray(a.obsm["X_pca"])
rsc.pp.neighbors(a, use_rep="X_pca")
rsc.tl.umap(a, random_state=0)
E.obsm["X_umap_unintegrated"] = np.asarray(a.obsm["X_umap"])
log("unintegrated done")
rsc.pp.harmony_integrate(a, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony")
E.obsm["X_pca_harmony"] = np.asarray(a.obsm["X_pca_harmony"])
rsc.pp.neighbors(a, use_rep="X_pca_harmony")
rsc.tl.umap(a, random_state=0)
rsc.tl.leiden(a, resolution=1.0, random_state=0)
E.obsm["X_umap_harmony"] = np.asarray(a.obsm["X_umap"])
E.obs["leiden_harmony_r1"] = a.obs["leiden"].astype(str).values
log("harmony done")
E.write_h5ad(f"{D}/embeddings_harmony.h5ad")
log("written")
