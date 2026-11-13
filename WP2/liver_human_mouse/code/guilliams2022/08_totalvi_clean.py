"""Step 08 - TotalVI on the cleaned counts (07_clean, 321,059 cells/nuclei), same settings as the draft (03_integrate.py):
4000 HVG (seurat flavour, batch-aware by sample) on log-normalised cleaned counts, TOTALVI n_latent 20, encode_covariates,
lr 1e-3, 400 epochs, train 0.9, batch = sample, batch_size 1024 + 8 loader workers (draft deviation kept for comparability).
Counts = X of the clean file (decontX (1,300) for 13 in vivo samples, raw otherwise; integers). Proteins = obsm adt
(all-zero in non-CITE batches -> masked by scvi-tools).
Then on GPU: neighbours, UMAP, Leiden r1/r2/r3; UMAP panel of the curated draft labels.
Output <OUT>/08_totalvi_clean/: embeddings.h5ad (obs + X_totalvi, X_umap_totalvi, leiden), totalvi_model/, history, figures/."""
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

OUT = "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/08_totalvi_clean"
os.makedirs(f"{D}/figures", exist_ok=True)
t0 = time.time()


def log(m):
    print(f"[{(time.time() - t0) / 60:6.1f} min] {m}", flush=True)


a = ad.read_h5ad(f"{OUT}/07_clean/guilliams2022_fig1b_clean.h5ad")
a.layers.pop("counts")  # raw; not used here
a.layers["clean_counts"] = a.X.copy()
log(f"loaded {a.shape}")
b = a.copy()
sc.pp.normalize_total(b, target_sum=1e4)
sc.pp.log1p(b)
sc.pp.highly_variable_genes(b, n_top_genes=4000, flavor="seurat", batch_key="sample")
hvg = b.var["highly_variable"].values
del b
log(f"HVG {hvg.sum()}")
pd.Series(a.var_names[hvg]).to_csv(f"{D}/hvg.tsv", sep="\t", index=False)

import scvi  # noqa: E402
import torch  # noqa: E402
from lightning.pytorch.callbacks import Callback  # noqa: E402

scvi.settings.seed = 0
t = a[:, hvg].copy()
t.obsm["protein_expression"] = np.asarray(a.obsm["adt"])
scvi.model.TOTALVI.setup_anndata(t, layer="clean_counts", batch_key="sample", protein_expression_obsm_key="protein_expression")
m = scvi.model.TOTALVI(t, n_latent=20, encode_covariates=True)
log(f"TOTALVI set up; cuda={torch.cuda.is_available()}")


class EpochLog(Callback):
    def __init__(self):
        self.t = time.time()

    def on_train_epoch_end(self, trainer, pl_module):
        e = trainer.current_epoch + 1
        if e % 5 == 0 or e == 1:
            dt = time.time() - self.t
            print(f"epoch {e}/{trainer.max_epochs} {dt / e:.1f} s/epoch, eta {dt / e * (trainer.max_epochs - e) / 60:.0f} min", flush=True)


scvi.settings.dl_num_workers = 8
m.train(max_epochs=400, lr=1e-3, batch_size=1024, train_size=0.9, accelerator="gpu", devices=1,
        enable_progress_bar=False, early_stopping=False, callbacks=[EpochLog()])
log("TOTALVI trained")
m.save(f"{D}/totalvi_model", overwrite=True)
pd.concat([v for v in m.history.values() if isinstance(v, pd.DataFrame)], axis=1).to_csv(f"{D}/totalvi_history.tsv", sep="\t")
Z = m.get_latent_representation()
del m, t
torch.cuda.empty_cache()
E = ad.AnnData(obs=a.obs.copy(), obsm={"X_totalvi": Z.astype(np.float32)})
for k in E.obs.columns:
    if E.obs[k].dtype == object:
        E.obs[k] = E.obs[k].astype(str)
E.write_h5ad(f"{D}/embeddings.h5ad")
log("latent written; UMAP/Leiden/panel in 08p_post_gpu.py")
