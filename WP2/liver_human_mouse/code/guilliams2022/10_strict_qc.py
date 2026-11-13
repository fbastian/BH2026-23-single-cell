"""Step 10 - strict QC approximating the authors' (user go 2026-11-12). Input: 07_clean (321,059) + 08 clean TotalVI latent.
Removed (all reasons recorded in obs strict_remove_reason, first reason wins in this order):
  1. scDblFinder_doublet  per-cell scDblFinder call (deviation from the MASLD consensus rule, to match the authors' strict
                          doublet removal: "eliminating these potential doublets", Cell 2022 p.382)
  2. dx_lt200             < 200 UMI after decontX (decontX samples, neutrophils exempt)
  3. bridge               > 50 % of the 15 nearest neighbours (clean TotalVI latent) carry another annotation_curated label
                          (intermediate profile between cell types; 'unresolved' neighbours count as different)
  4. sort_mismatch        cell type not expected for the sample's FACS gate (GEO sample titles, GSE192740 SOFT):
                          CD45+ -> immune + HsPCs; CD45- -> hepatocyte, endothelial, fibroblast, cholangiocyte (checked 2026-11-12: removed cells in CD45+ gates are Ptprc RNA 0-5 % positive vs 78-98 % in kept immune cells);
                          myeloid / monos+macs / DCs -> myeloid; hepatocytes -> hepatocyte; total live (CS144) and nuclei -> all.
                          'unresolved' cells are never removed for sort mismatch.
Writes <OUT>/10_strict/: guilliams2022_fig1b_strict.h5ad (kept cells; X = clean counts, layers counts = raw), cells_strict.parquet
(all 321,059 with reasons), per_sample_strict.tsv, removal_by_type.tsv, sample_sort_gate.tsv."""
import gzip
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
from cuml.neighbors import NearestNeighbors

OUT = "/mnt/masld_vault/cross_species/guilliams2022"
D = f"{OUT}/10_strict"
os.makedirs(D, exist_ok=True)

# sort gate per sample from GEO titles
title = {}
with gzip.open("/mnt/masld_vault/geo/GSE192740/_meta/family.soft.gz", "rt") as fh:
    gsm = None
    for line in fh:
        if line.startswith("^SAMPLE"):
            gsm = line.split("=")[1].strip()
        elif line.startswith("!Sample_title") and gsm:
            title[gsm] = line.split("=", 1)[1].strip()
st = pd.read_csv(f"{OUT}/01_raw/sample_table.tsv", sep="\t")
st["title"] = st.gsm.map(title)


def gate(t):
    t = t.lower()
    if "nuclei" in t:
        return "nuclei"
    if "hepatocytes" in t:
        return "hepatocytes"
    if re.search(r"cd45\+", t):
        return "CD45+"
    if re.search(r"cd45-", t):
        return "CD45-"
    if "myeloid" in t or "monos+macs" in t or re.search(r"(^|[^a-z])dcs([^a-z]|$)", t):
        return "myeloid"
    return "total_live"


st["gate"] = st.title.map(gate)
st[["sample", "gsm", "gate", "title"]].to_csv(f"{D}/sample_sort_gate.tsv", sep="\t", index=False)
print(st.gate.value_counts().to_string(), flush=True)
print(st[["sample", "gate", "title"]].to_string(index=False), flush=True)

IMMUNE = {"KC", "Monocyte", "T cell", "pDCs", "B cell", "cDC2", "cDC1", "ILC1", "Neutrophil", "NK cell", "Mig. cDCs",
          "Basophil", "HsPCs"}
MYELOID = {"KC", "Monocyte", "pDCs", "cDC2", "cDC1", "Neutrophil", "Mig. cDCs", "Basophil"}
NONIMM = {"Hepatocyte", "Endothelial", "Fibroblast", "Cholangiocyte"}
ALLOWED = {"CD45+": IMMUNE, "CD45-": NONIMM, "myeloid": MYELOID, "hepatocytes": {"Hepatocyte"}}

a = ad.read_h5ad(f"{OUT}/07_clean/guilliams2022_fig1b_clean.h5ad")
E = ad.read_h5ad(f"{OUT}/08_totalvi_clean/embeddings.h5ad")
assert (E.obs_names == a.obs_names).all()
o = a.obs
o["sort_gate"] = o["sample"].astype(str).map(st.set_index("sample").gate).values
lab = o["annotation_curated"].astype(str).values
Z = E.obsm["X_totalvi"].astype(np.float32)
_, idx = NearestNeighbors(n_neighbors=16).fit(Z).kneighbors(Z)
idx = np.asarray(idx)[:, 1:]
o["bridge_type"] = (lab[idx] != lab[:, None]).mean(1)
mism = np.array([g in ALLOWED and l != "unresolved" and l not in ALLOWED[g] for g, l in zip(o.sort_gate.values, lab)])
reason = np.full(len(o), "", dtype=object)
for name, m in [("scDblFinder_doublet", o.scDblFinder_doublet.values.astype(bool)), ("dx_lt200", o.dx_lt200.values.astype(bool)),
                ("bridge", o.bridge_type.values > 0.5), ("sort_mismatch", mism)]:
    reason[(reason == "") & m] = name
o["strict_remove_reason"] = pd.Categorical(np.where(reason == "", "kept", reason))
o["strict_keep"] = reason == ""
o[["sample", "sort_gate", "suspension", "annotation_curated", "bridge_type", "strict_remove_reason", "strict_keep"]].to_parquet(
    f"{D}/cells_strict.parquet")
print(o.strict_remove_reason.value_counts().to_string(), "\nkept", int(o.strict_keep.sum()), flush=True)
ps = pd.crosstab(o["sample"], o.strict_remove_reason)
ps.insert(0, "gate", st.set_index("sample").gate.reindex(ps.index).values)
ps.to_csv(f"{D}/per_sample_strict.tsv", sep="\t")
bt = pd.crosstab(o.annotation_curated, o.strict_remove_reason)
bt.to_csv(f"{D}/removal_by_type.tsv", sep="\t")
print(bt.to_string(), flush=True)
k = a[o.strict_keep.values].copy()
for c in k.obs.columns:
    if k.obs[c].dtype == object:
        k.obs[c] = k.obs[c].astype(str)
k.uns["note"] = ("STRICT 2026-11-12: 07_clean minus per-cell scDblFinder doublets, <200 UMI after decontX, bridge cells (>50 % of 15 NN "
                 "with another curated label), sort-gate mismatches (GEO titles). See 10_strict/cells_strict.parquet.")
k.write_h5ad(f"{D}/guilliams2022_fig1b_strict.h5ad", compression="gzip")
print("written", k.shape, flush=True)
