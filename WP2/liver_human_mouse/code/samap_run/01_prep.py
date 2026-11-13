"""SAMap prep (2026-11-12): HLiCA v1.0 x Guilliams 2022 mouse atlas (annotation of record, 13b).
Subsample stratified by cell type (seed 0): mouse cap 4,000 per annotation_aor ('unresolved' excluded); human cap chosen so the
total is ~50k over the 48 author_cell_type labels. X = counts (mouse: clean counts = decontX (1,300) for 13 in vivo samples, raw
otherwise; human: raw.X UMI). var_names -> SAMap names (refs/mappings/{hs,mm}.gene_names.tsv, Symbol or Symbol__EnsemblID);
genes absent from the tables dropped. Outputs /mnt/masld_vault/cross_species/samap_run/input/{hs,mm}_samap_ready.h5ad + subsample tables."""
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
MAPS = "/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/refs/mappings"
O = "/mnt/masld_vault/cross_species/samap_run/input"
rng = np.random.default_rng(0)
def strat(lab, cap):
    idx = []
    for t, ii in pd.Series(np.arange(len(lab))).groupby(lab.values).groups.items():
        ii = np.asarray(ii); idx.append(ii if len(ii) <= cap else rng.choice(ii, cap, replace=False))
    return np.sort(np.concatenate(idx))
def rename(a, tab):
    g = pd.read_csv(tab, sep="\t").set_index("ensembl_gene_id").samap_gene
    keep = a.var_names.isin(g.index); a = a[:, keep].copy(); a.var["ensembl_gene_id"] = a.var_names.astype(str)
    a.var_names = g.loc[a.var_names].astype(str).values; assert a.var_names.is_unique; a.var.index.name = "samap_gene"
    print("genes kept", keep.sum(), "of", len(keep), flush=True); return a
# mouse
m = ad.read_h5ad("/mnt/masld_vault/cross_species/guilliams2022/13_annotation_of_record/guilliams2022_fig1b_aor.h5ad")
m = m[m.obs.annotation_aor.astype(str) != "unresolved"]
ix = strat(m.obs.annotation_aor.astype(str), 4000)
mm = ad.AnnData(X=sp.csr_matrix(m.X[ix]).astype(np.float32), obs=m.obs.iloc[ix][["sample", "suspension", "digestion", "sort_gate", "annotation_aor", "aor_flag"]].copy(),
                var=m.var[["symbol"]].copy())
mm.obs["annotation_aor"] = mm.obs.annotation_aor.astype(str).astype("category"); mm.obsm["X_umap_mouse"] = m.obsm["X_umap"][ix]
del m
mm = rename(mm, f"{MAPS}/mm.gene_names.tsv")
print(mm.obs.annotation_aor.value_counts().to_string(), "\nmouse", mm.shape, flush=True)
mm.obs.annotation_aor.value_counts().to_csv(f"{O}/mm_subsample_counts.tsv", sep="\t")
mm.write_h5ad(f"{O}/mm_samap_ready.h5ad")
# human
h = ad.read_h5ad("/mnt/masld_vault/geo/cxg-059202e1-1f1b-483f-9151-f3a25a380c39/HLiCA_human_liver_all_cells.h5ad", backed="r")
lab = h.obs.author_cell_type.astype(str); vc = lab.value_counts()
cap = next(c for c in range(200, 20000, 50) if np.minimum(vc.values, c).sum() >= 50000)
print("human cap per type", cap, flush=True)
ix = strat(lab, cap)
assert (h.raw.var_names == h.var_names).all()
X = h.raw.X[ix] if not hasattr(h.raw.X, "shape") or True else None
X = sp.csr_matrix(X).astype(np.float32)
hs = ad.AnnData(X=X, obs=h.obs.iloc[ix][["author_cell_type", "cell_type", "suspension_type", "STUDY", "donor_id", "assay", "condition"]].copy(),
                var=h.var[["feature_name", "feature_biotype"]].copy())
hs.obsm["X_umap_hlica"] = np.asarray(h.obsm["X_umap"][ix])
for c in hs.obs.columns:
    if hs.obs[c].dtype.name == "category": hs.obs[c] = hs.obs[c].cat.remove_unused_categories()
hs = rename(hs, f"{MAPS}/hs.gene_names.tsv")
print(hs.obs.author_cell_type.value_counts().to_string(), "\nhuman", hs.shape, flush=True)
hs.obs.author_cell_type.value_counts().to_csv(f"{O}/hs_subsample_counts.tsv", sep="\t")
hs.write_h5ad(f"{O}/hs_samap_ready.h5ad")
print("done", flush=True)
