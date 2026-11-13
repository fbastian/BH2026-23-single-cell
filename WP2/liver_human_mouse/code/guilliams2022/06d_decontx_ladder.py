"""Step 06d - decontX strength ladder on one sample (MASLD decontx_titration idea): default (estimated delta) vs fixed priors
(1,100), (1,300), (1,1000), estimateDelta FALSE; z = TotalVI Leiden r1; varGenes 5000; seed 12345.
Metrics: Alb removed from non-epithelial draft types, Alb kept in hepatocytes, UMIs kept, cells < 200 UMI.
Usage: 06d_decontx_ladder.py <SPLITDIR> <RESDIR> <sample>  -> <RESDIR>/<sample>.ladder.tsv"""
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
alb = np.flatnonzero(a.var["symbol"].astype(str).values[gk] == "Alb")[0]
typ = a.obs["annotation_draft"].astype(str).values
hep, nonepi = typ == "Hepatocyte", ~np.isin(typ, ["Hepatocyte", "Cholangiocyte"])
a0 = np.asarray(X[:, gk[alb]].todense()).ravel(); u0 = np.asarray(X.sum(1)).ravel()
rows = []
for name, arg in [("default_est", ""), ("1_100", ", delta = c(1, 100), estimateDelta = FALSE"),
                  ("1_300", ", delta = c(1, 300), estimateDelta = FALSE"), ("1_1000", ", delta = c(1, 1000), estimateDelta = FALSE")]:
    ro.r(f"res <- decontX(m, z = as.character(zl), varGenes = 5000, seed = 12345, verbose = FALSE{arg}); "
         f"dc <- res$decontXcounts; cont <- res$contamination; ab <- as.numeric(dc[{alb + 1}, ]); us <- colSums(dc)")
    a1 = np.asarray(ro.r("ab")); u1 = np.asarray(ro.r("us")); cont = np.asarray(ro.r("cont"))
    rows.append(dict(sample=S, setting=name, n=len(a1), median_contamination=np.median(cont),
                     nonepi_Alb_cp10k_raw=(a0[nonepi] / np.maximum(u0[nonepi], 1) * 1e4).mean(),
                     nonepi_Alb_cp10k_dx=(a1[nonepi] / np.maximum(u1[nonepi], 1) * 1e4).mean(),
                     nonepi_Alb_removed=1 - a1[nonepi].sum() / max(a0[nonepi].sum(), 1),
                     hep_Alb_kept=a1[hep].sum() / a0[hep].sum() if a0[hep].sum() else np.nan,
                     hep_umis_kept=u1[hep].sum() / u0[hep].sum() if hep.any() else np.nan,
                     nonepi_umis_kept=u1[nonepi].sum() / u0[nonepi].sum(), lt200=int((u1 < 200).sum())))
    print(rows[-1], flush=True)
os.makedirs(RES, exist_ok=True)
pd.DataFrame(rows).to_csv(f"{RES}/{S}.ladder.tsv", sep="\t", index=False, float_format="%.4g")
