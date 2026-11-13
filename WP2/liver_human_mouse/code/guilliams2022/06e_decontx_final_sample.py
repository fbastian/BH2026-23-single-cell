"""Step 06e - final decontX for one tier-1 sample: fixed prior delta = (1, 300), estimateDelta FALSE, z = TotalVI Leiden r1,
varGenes 5000, seed 12345 (chosen from the 06d ladder, user go 2026-11-11). Draft-QC kept cells only.
Usage: 06e_decontx_final_sample.py <SPLITDIR> <RESDIR> <sample>  -> <RESDIR>/<sample>.dx300.npz, <sample>.dx300.parquet"""
import os
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.simplefilter("ignore")
SPLIT, RES, S = sys.argv[1:4]
os.environ.setdefault("R_HOME", f"{sys.prefix}/lib/R")
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import numpy2ri  # noqa: E402

ro.r("suppressPackageStartupMessages({library(celda); library(Matrix)})")
CONV = ro.default_converter + numpy2ri.converter
a = ad.read_h5ad(f"{SPLIT}/{S}.h5ad")
a = a[a.obs["keep"].values].copy()
X = sp.csr_matrix(a.X)
gk = np.flatnonzero(np.asarray((X > 0).sum(0)).ravel() > 0)
m = sp.csc_matrix(X[:, gk].T); m.sort_indices()
with CONV.context():
    ro.globalenv["._i"] = m.indices.astype(np.int32); ro.globalenv["._p"] = m.indptr.astype(np.int32)
    ro.globalenv["._x"] = m.data.astype(np.float64); ro.globalenv["._d"] = np.asarray(m.shape, dtype=np.int32)
    ro.globalenv["zl"] = a.obs["z"].astype(str).values
ro.r('m <- new("dgCMatrix", i=as.integer(._i), p=as.integer(._p), x=as.numeric(._x), Dim=as.integer(._d))')
ro.r("res <- decontX(m, z = as.character(zl), delta = c(1, 300), estimateDelta = FALSE, varGenes = 5000, seed = 12345, "
     "verbose = FALSE); dc <- drop0(round(res$decontXcounts)); cont <- res$contamination; rm(res)")
with CONV.context():
    i = np.asarray(ro.r("dc@i")); p = np.asarray(ro.r("dc@p")); x = np.asarray(ro.r("dc@x"))
cg = sp.csc_matrix((x.astype(np.float32), i, p), shape=(len(gk), a.n_obs)).T.tocsr()
corr = sp.csr_matrix((cg.data, gk[cg.indices], cg.indptr), shape=X.shape, dtype=np.float32)
sp.save_npz(f"{RES}/{S}.dx300.npz", corr)
pd.DataFrame({"decontX_contamination_final": np.asarray(ro.r("cont"))}, index=a.obs_names).to_parquet(f"{RES}/{S}.dx300.parquet")
print(S, a.n_obs, f"UMIs kept {corr.sum() / X.sum():.3f}", flush=True)
