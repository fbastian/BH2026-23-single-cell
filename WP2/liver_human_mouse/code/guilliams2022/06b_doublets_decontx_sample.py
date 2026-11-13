"""Step 06b - one sample: scDblFinder on all cells (MASLD doublets_v6 settings: defaults, raw counts, seed 12345), then
decontX on the draft-QC kept cells (MASLD fq_08 settings: z = clusters, varGenes 5000, seed 12345; here z = TotalVI Leiden r1
and decontX defaults: delta c(10,10), estimateDelta TRUE - no ambient tiers exist for this dataset yet).
sc-preprocessing env (R scDblFinder + celda via rpy2).  Usage: 06b_doublets_decontx_sample.py <SPLITDIR> <RESDIR> <sample>
Writes <RESDIR>/<sample>.cells.parquet (scDblFinder score/class, decontX contamination) and <sample>.decontx.npz (kept cells)."""
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.simplefilter("ignore")
SPLIT, RES, S = sys.argv[1:4]
SEED = 12345
os.environ.setdefault("R_HOME", f"{sys.prefix}/lib/R")
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import numpy2ri  # noqa: E402

ro.r("suppressPackageStartupMessages({library(scDblFinder); library(SingleCellExperiment); library(celda); library(Matrix)})")
CONV = ro.default_converter + numpy2ri.converter


def push(m, name):
    m = sp.csc_matrix(m)
    m.sort_indices()
    with CONV.context():
        ro.globalenv["._i"] = m.indices.astype(np.int32)
        ro.globalenv["._p"] = m.indptr.astype(np.int32)
        ro.globalenv["._x"] = m.data.astype(np.float64)
        ro.globalenv["._d"] = np.asarray(m.shape, dtype=np.int32)
    ro.r(f'{name} <- new("dgCMatrix", i=as.integer(._i), p=as.integer(._p), x=as.numeric(._x), Dim=as.integer(._d))')


def pull(name, shape):
    with CONV.context():
        i = np.asarray(ro.r(f"{name}@i")); p = np.asarray(ro.r(f"{name}@p")); x = np.asarray(ro.r(f"{name}@x"))
    return sp.csc_matrix((x.astype(np.float32), i, p), shape=shape)


t0 = time.time()
a = ad.read_h5ad(f"{SPLIT}/{S}.h5ad")
X = sp.csr_matrix(a.X)
assert np.array_equal(X.data, np.round(X.data)), "not integer counts"
res = pd.DataFrame(index=a.obs_names, data={"scDblFinder_score": np.nan, "scDblFinder_class": "not_run",
                                            "decontX_contamination": np.nan})
push(X.T, "m")
ro.r(f"set.seed({SEED}); sce <- scDblFinder(SingleCellExperiment(list(counts = m)), verbose = FALSE); "
     "dsc <- as.numeric(sce$scDblFinder.score); dcl <- as.character(sce$scDblFinder.class); rm(sce, m); invisible(gc())")
res["scDblFinder_score"] = np.asarray(ro.r("dsc"))
res["scDblFinder_class"] = list(ro.r("dcl"))
t1 = time.time()

k = np.flatnonzero(a.obs["keep"].values)
sub = X[k]
gk = np.flatnonzero(np.asarray((sub > 0).sum(0)).ravel() > 0)
push(sub[:, gk].T, "m")
z = a.obs["z"].values[k].astype(str)
with CONV.context():
    ro.globalenv["zl"] = np.asarray(z, dtype=str)
zarg = "as.character(zl)" if len(set(z)) > 1 else "NULL"
ro.r(f"res <- decontX(m, z = {zarg}, varGenes = 5000, seed = {SEED}, verbose = FALSE); "
     "dc <- drop0(round(res$decontXcounts)); cont <- res$contamination; est <- res$estimates; rm(res)")
cg = pull("dc", (len(gk), len(k))).T.tocsr()
corr = sp.csr_matrix((cg.data, gk[cg.indices], cg.indptr), shape=sub.shape, dtype=np.float32)
res.iloc[k, 2] = np.asarray(ro.r("cont"))
os.makedirs(RES, exist_ok=True)
sp.save_npz(f"{RES}/{S}.decontx.npz", corr)
res.to_parquet(f"{RES}/{S}.cells.parquet")
print(f"{S}: {len(res)} cells, scDblFinder doublets {(res.scDblFinder_class == 'doublet').mean():.3f} ({t1 - t0:.0f}s); "
      f"decontX on {len(k)}: median contamination {np.nanmedian(res.decontX_contamination):.3f}, "
      f"UMIs kept {corr.sum() / sub.sum():.3f} ({time.time() - t1:.0f}s)", flush=True)
