"""Step 06a - split the raw Fig 1B object (all 389,056 cells) into one h5ad per sample on local disk for the
scDblFinder/decontX workers. obs carries keep (draft QC), scrublet_doublet and leiden_totalvi_r1 (kept cells only).
Usage: 06a_split_samples.py <OUT> <SPLITDIR>"""
import os
import sys

import anndata as ad
import pandas as pd

OUT, SPLIT = sys.argv[1], sys.argv[2]
os.makedirs(SPLIT, exist_ok=True)
a = ad.read_h5ad(f"{OUT}/01_raw/guilliams2022_fig1b_raw.h5ad")
f = pd.read_parquet(f"{OUT}/02_qc/guilliams2022_fig1b_qc_flags.parquet")
assert (f.index == a.obs_names).all()
E = pd.read_parquet(f"{OUT}/05_annotation_draft/cells_annotation_draft.parquet")
a.obs["keep"] = f["keep"].values
a.obs["scrublet_doublet"] = f["scrublet_doublet"].values
a.obs["scrublet_score"] = f["doublet_score"].values
a.obs["z"] = E["leiden_totalvi_r1"].reindex(a.obs_names).astype(str).values
a.obs["annotation_draft"] = E["annotation_draft"].reindex(a.obs_names).astype(str).values
assert (a.obs.loc[a.obs.keep, "z"] != "nan").all()
a.obsm.clear()
for s, idx in a.obs.groupby("sample", observed=True).indices.items():
    a[idx].copy().write_h5ad(f"{SPLIT}/{s}.h5ad")
    print(s, len(idx), flush=True)
