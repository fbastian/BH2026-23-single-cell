"""Step 02 - per-sample QC following the authors' 3a_SCT.R (scater isOutlier, 3 MADs), plus Scrublet doublets.

Authors (3a_SCT.R): per sample, outliers at 3 MADs on log(n genes) low+high, log(n UMI) low+high, % mito (^mt-) high.
They then removed doublet / contaminating / proliferating clusters by hand per sample; we stand in Scrublet
(per sample, scanpy defaults) for the doublet part and record everything as flags before removing.

Input : <OUT>/01_raw/guilliams2022_fig1b_raw.h5ad
Output: <OUT>/02_qc/guilliams2022_fig1b_qc_flags.parquet  (all cells, all flags)
        <OUT>/02_qc/guilliams2022_fig1b_qc.h5ad            (kept cells only, raw counts)
        <OUT>/02_qc/qc_per_sample.tsv
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
os.makedirs(f"{OUT}/02_qc", exist_ok=True)
NMAD = 3

a = ad.read_h5ad(f"{OUT}/01_raw/guilliams2022_fig1b_raw.h5ad")
mito = a.var["symbol"].str.startswith("mt-").values
a.obs["pct_mito"] = 100 * np.asarray(a.X[:, mito].sum(1)).ravel() / np.maximum(a.obs["n_counts"].values, 1)


def mad_out(x, lower, upper, log):
    x = np.log1p(x) if log else np.asarray(x, float)
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))  # scater isOutlier uses the scaled MAD
    lo = (x < med - NMAD * mad) if lower else np.zeros(len(x), bool)
    hi = (x > med + NMAD * mad) if upper else np.zeros(len(x), bool)
    return lo, hi


try:  # GPU Scrublet (rapids-singlecell), per sample via batch_key; CPU scanpy fallback below
    import cupy as cp
    import rapids_singlecell as rsc
    parts = []
    for s, idx in a.obs.groupby("sample", observed=True).indices.items():  # one sample at a time fits a T4
        g = ad.AnnData(X=a.X[idx].astype(np.float32), obs=a.obs.iloc[idx][["sample"]].copy())
        rsc.get.anndata_to_GPU(g)
        rsc.pp.scrublet(g, random_state=0)
        parts.append(g.obs[["doublet_score", "predicted_doublet"]].copy())
        del g
        cp.get_default_memory_pool().free_all_blocks()
    GPU_DBL = pd.concat(parts).reindex(a.obs_names)
    print("scrublet on GPU done", flush=True)
except Exception as ex:
    print(f"GPU scrublet unavailable ({ex}); CPU per sample", flush=True)
    GPU_DBL = None

flags = []
for s, idx in a.obs.groupby("sample", observed=True).indices.items():
    o = a.obs.iloc[idx]
    gl, gh = mad_out(o["n_genes"].values, True, True, True)
    ul, uh = mad_out(o["n_counts"].values, True, True, True)
    _, mh = mad_out(o["pct_mito"].values, False, True, False)
    f = pd.DataFrame({"genes_low": gl, "genes_high": gh, "umi_low": ul, "umi_high": uh, "mito_high": mh},
                     index=o.index)
    if GPU_DBL is not None:
        f["doublet_score"] = GPU_DBL["doublet_score"].iloc[idx].values
        f["scrublet_doublet"] = GPU_DBL["predicted_doublet"].iloc[idx].astype(bool).values
        flags.append(f)
        print(f"{s:8s} n={len(o):6d} mad_out={int((f.iloc[:, :5].any(axis=1)).sum()):5d} "
              f"dbl={int(f['scrublet_doublet'].sum()):5d}", flush=True)
        continue
    sub = a[idx].copy()
    try:
        sc.pp.scrublet(sub, random_state=0)
        f["doublet_score"] = sub.obs["doublet_score"].values
        f["scrublet_doublet"] = sub.obs["predicted_doublet"].astype(bool).values
    except Exception as ex:  # tiny samples can fail
        print(f"scrublet failed for {s}: {ex}", flush=True)
        f["doublet_score"] = np.nan
        f["scrublet_doublet"] = False
    flags.append(f)
    print(f"{s:8s} n={len(o):6d} mad_out={int((f.iloc[:, :5].any(axis=1)).sum()):5d} "
          f"dbl={int(f['scrublet_doublet'].sum()):5d}", flush=True)

F = pd.concat(flags).reindex(a.obs_names)
F["mad_outlier"] = F[["genes_low", "genes_high", "umi_low", "umi_high", "mito_high"]].any(axis=1)
F["keep"] = ~(F["mad_outlier"] | F["scrublet_doublet"])
F = F.join(a.obs[["sample", "digestion", "suspension", "cite_seq", "n_counts", "n_genes", "pct_mito"]])
F.to_parquet(f"{OUT}/02_qc/guilliams2022_fig1b_qc_flags.parquet")

per = F.groupby("sample", observed=True).agg(digestion=("digestion", "first"), n=("keep", "size"),
                                             mad_outlier=("mad_outlier", "sum"),
                                             doublet=("scrublet_doublet", "sum"), kept=("keep", "sum"))
per.to_csv(f"{OUT}/02_qc/qc_per_sample.tsv", sep="\t")
print(per.to_string())
print("by digestion:\n", F.groupby("digestion", observed=True)["keep"].agg(["size", "sum"]).to_string())
print(f"TOTAL kept {int(F['keep'].sum())} of {len(F)} (paper after QC: 185,894)", flush=True)

for c in ["pct_mito"]:
    a.obs[c] = a.obs[c].astype(np.float32)
a.obs = a.obs.join(F[["doublet_score"]])
k = a[F["keep"].values].copy()
k.write_h5ad(f"{OUT}/02_qc/guilliams2022_fig1b_qc.h5ad", compression="gzip")
print("written", k.shape, flush=True)
