"""Step 03p - neighbours, UMAP and Leiden (r=1) on the TotalVI latent, on GPU (rapids-singlecell).
Updates <OUT>/03_integration/embeddings.h5ad with X_umap_totalvi and leiden_totalvi_r1."""
import sys

import anndata as ad
import numpy as np
import rapids_singlecell as rsc

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
p = f"{OUT}/03_integration/embeddings.h5ad"
E = ad.read_h5ad(p)
t = ad.AnnData(obs=E.obs[[]].copy(), obsm={"X_totalvi": E.obsm["X_totalvi"].astype(np.float32)})
rsc.pp.neighbors(t, use_rep="X_totalvi")
rsc.tl.umap(t, random_state=0)
rsc.tl.leiden(t, resolution=1.0, random_state=0)
E.obsm["X_umap_totalvi"] = np.asarray(t.obsm["X_umap"])
E.obs["leiden_totalvi_r1"] = t.obs["leiden"].astype(str).values
E.write_h5ad(p)
print("written", E.shape, E.obs["leiden_totalvi_r1"].nunique(), "clusters", flush=True)
