"""Step 07 - apply the manual cluster overrides (07_clean/nuclei_overrides.tsv) to the three clean files.
Adds obs annotation_curated (= annotation_clean unless overridden) and annotation_flag; rewrites the files in place (tmp + rename)."""
import os

import anndata as ad
import pandas as pd

D = "/mnt/masld_vault/cross_species/guilliams2022/07_clean"
ov = pd.read_csv(f"{D}/nuclei_overrides.tsv", sep="\t", comment="#", dtype=str).fillna("")
ov = ov.set_index("leiden_totalvi_r3")
rows = []
for f in ("guilliams2022_fig1b_clean.h5ad", "guilliams2022_fig1b_clean_cells.h5ad", "guilliams2022_fig1b_clean_nuclei.h5ad"):
    a = ad.read_h5ad(f"{D}/{f}")
    cl = a.obs["leiden_totalvi_r3"].astype(str)
    cur = a.obs["annotation_clean"].astype(str).copy()
    flag = pd.Series("", index=a.obs_names)
    hit = cl.isin(ov.index)
    cur[hit] = ov.loc[cl[hit], "annotation_curated"].values
    flag[hit] = ov.loc[cl[hit], "flag"].values
    a.obs["annotation_curated"] = pd.Categorical(cur)
    a.obs["annotation_flag"] = pd.Categorical(flag)
    a.uns["annotation_overrides"] = "07_clean/nuclei_overrides.tsv (user decision 2026-11-11), applied per TotalVI draft Leiden r3 cluster"
    a.write_h5ad(f"{D}/{f}.tmp", compression="gzip")
    os.replace(f"{D}/{f}.tmp", f"{D}/{f}")
    t = a.obs.groupby(["annotation_curated", "suspension"], observed=True).size().unstack(fill_value=0)
    print(f, a.shape, "overridden", int(hit.sum()), flush=True)
    if f.endswith("clean.h5ad"):
        t["total"] = t.sum(1)
        t.sort_values("total", ascending=False).to_csv(f"{D}/composition_curated.tsv", sep="\t")
        print(t.sort_values("total", ascending=False).to_string(), flush=True)
