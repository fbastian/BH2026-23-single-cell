"""Step 03 - integration of the Fig 1B cells: unintegrated PCA, Harmony baseline, TotalVI (authors' method).

Authors (4c_TotalVI_CiteSeq+non-CiteSeq.ipynb, scvi 0.6 API): union of genes/proteins, 4000 HVG
(seurat_v2, batch-corrected), proteins of non-CITE samples set to 0 and masked per batch, TOTALVI n_latent=20,
encoder_batch=True, lr=1e-3, 400 epochs, train 90 %. Batch = sample.
Here (scvi-tools 1.5): same settings; scvi-tools masks proteins that are all-zero within a batch automatically.

Usage: 03_integrate.py <OUT> [harmony|totalvi|both]
Output: <OUT>/03_integration/
  embeddings.h5ad  obs + obsm X_pca, X_pca_harmony, X_totalvi, X_umap_<method>; no expression matrix
  totalvi_model/   saved model
"""
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
MODE = sys.argv[2] if len(sys.argv) > 2 else "both"
D = f"{OUT}/03_integration"
os.makedirs(D, exist_ok=True)
t0 = time.time()


def log(m):
    print(f"[{(time.time() - t0) / 60:6.1f} min] {m}", flush=True)


a = ad.read_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad")
a.layers["counts"] = a.X.copy()
log(f"loaded {a.shape}")

# HVG: 4000, batch-aware, Seurat flavour on log-normalised data (authors: seurat_v2 with batch correction)
b = a.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)
sc.pp.highly_variable_genes(b, n_top_genes=4000, flavor="seurat", batch_key="sample")
hvg = b.var["highly_variable"].values
log(f"HVG {hvg.sum()}")

emb_path = f"{D}/embeddings.h5ad"
E = ad.read_h5ad(emb_path) if os.path.exists(emb_path) else ad.AnnData(obs=a.obs.copy())
assert list(E.obs_names) == list(a.obs_names)

if MODE in ("harmony", "both"):
    c = b[:, hvg].copy()
    sc.pp.scale(c, max_value=10)
    sc.pp.pca(c, n_comps=50, random_state=0)
    E.obsm["X_pca"] = c.obsm["X_pca"]
    sc.pp.neighbors(c, use_rep="X_pca", random_state=0)
    sc.tl.umap(c, random_state=0)
    E.obsm["X_umap_unintegrated"] = c.obsm["X_umap"]
    log("unintegrated PCA/UMAP done")
    import harmonypy
    ho = harmonypy.run_harmony(c.obsm["X_pca"], c.obs, ["sample"], random_state=0)
    Z = np.asarray(ho.Z_corr)
    Z = Z.T if Z.shape[0] != c.n_obs else Z
    c.obsm["X_pca_harmony"] = Z
    E.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(c, use_rep="X_pca_harmony", random_state=0)
    sc.tl.umap(c, random_state=0)
    E.obsm["X_umap_harmony"] = c.obsm["X_umap"]
    log("harmony done")
    E.write_h5ad(emb_path)
    del c

if MODE in ("totalvi", "both"):
    import scvi
    import torch
    scvi.settings.seed = 0
    t = a[:, hvg].copy()
    t.obsm["protein_expression"] = a.obsm["adt"].copy()
    scvi.model.TOTALVI.setup_anndata(t, layer="counts", batch_key="sample",
                                     protein_expression_obsm_key="protein_expression")
    m = scvi.model.TOTALVI(t, n_latent=20, encode_covariates=True)
    log(f"TOTALVI set up; cuda={torch.cuda.is_available()}")
    from lightning.pytorch.callbacks import Callback

    class EpochLog(Callback):  # one line per 5 epochs with s/epoch and ETA
        def __init__(self):
            self.t = time.time()

        def on_train_epoch_end(self, trainer, pl_module):
            e = trainer.current_epoch + 1
            if e % 5 == 0 or e == 1:
                dt = time.time() - self.t
                print(f"epoch {e}/{trainer.max_epochs} {dt / e:.1f} s/epoch, "
                      f"eta {dt / e * (trainer.max_epochs - e) / 60:.0f} min", flush=True)

    # 2026-11-11 restart (user option 2): batch 1024 + 8 loader workers instead of scvi default 256 / 0 workers;
    # the 256 run was loader-bound (GPU ~38 %). Deviation from the authors' training settings, draft only.
    scvi.settings.dl_num_workers = 8
    m.train(max_epochs=400, lr=1e-3, batch_size=1024, train_size=0.9, accelerator="gpu", devices=1,
            enable_progress_bar=False, early_stopping=False, callbacks=[EpochLog()])
    log("TOTALVI trained")
    m.save(f"{D}/totalvi_model", overwrite=True)
    hist = pd.concat([v for v in m.history.values() if isinstance(v, pd.DataFrame)], axis=1)
    hist.to_csv(f"{D}/totalvi_history.tsv", sep="\t")
    E.obsm["X_totalvi"] = m.get_latent_representation()
    E.write_h5ad(emb_path)  # UMAP + Leiden on GPU in 03p_post_gpu.py
    log("latent written")

log("finished")
