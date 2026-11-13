"""Step 13b - annotation of record (AoR) for the strict object, decisions from 13a evidence (2026-11-12 CET).
Rule: label = Table S1 best type per Leiden r3 cluster (11_totalvi_strict, z within cell / nucleus clusters), except:
  - transitional DC: DC-block subcluster 3 (13a dc_sub; Siglech+ Ly6c2+ Cd209a+ Klk1- Ccr2+; profile 12b)
  - r3 12 (basophils + Siglech-low Ccr9+ pDC-like cells): per-cell curated label kept (pDCs / Basophil), flag mixed_r3_12
  - r3 80 (Dcn+ Krt19+ Cd34+): Mesothelial if Msln and Upk3b > 1 (not a Fig 1B type), else Fibroblast
  - unresolved (excluded from SAMap): r3 1 (hepatocyte+LSEC nuclei, scDblFinder 0.33), 61 (Alb+Stab2+Kdr, 442 UMI),
    84 (Alb+Dcn, 23 cells, one sample), 0 and 35 (low-UMI nuclei, weak Table S1 scores)
  - r3 29 -> cDC2 (S1 margin 0 vs Monocyte; H2-Ab1 4.8, Cx3cr1/Ccr2+) flag moDC_like; r3 65 -> ILC1 (S1); r3 82, 60 -> Hepatocyte
  - flags only: 39 KC_with_LSEC_signal, 63 B_with_myeloid_signal_high_UMI, 67 KC_with_RBC, 79 KC_low_UMI, 83 proliferating_single_sample
Writes <OUT>/13_annotation_of_record/: guilliams2022_fig1b_aor.h5ad (+ _cells, _nuclei; X clean counts, layers counts raw,
obsm X_totalvi + X_umap), cells_aor.tsv.gz, r3_labels_aor.tsv, composition_aor.tsv, figures/."""
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = "/mnt/masld_vault/cross_species/guilliams2022"; D = f"{OUT}/13_annotation_of_record"
E = ad.read_h5ad(f"{OUT}/11_totalvi_strict/embeddings.h5ad")
a = ad.read_h5ad(f"{OUT}/10_strict/guilliams2022_fig1b_strict.h5ad")
assert (a.obs_names == E.obs_names).all()
ev = pd.read_csv(f"{D}/r3_evidence.tsv", sep="\t", index_col=0); ev.index = ev.index.astype(str)
sub = pd.read_csv(f"{D}/dc_sub.tsv.gz", sep="\t", index_col=0).dc_sub.astype(str)
r3 = E.obs.leiden_clean_r3.astype(str)
n = a.copy(); sc.pp.normalize_total(n, target_sum=1e4); sc.pp.log1p(n)
n.var_names = n.var.symbol.astype(str).values; n.var_names_make_unique()
m80 = (r3 == "80").values
msl = {g: float(n[m80, g].X.mean()) for g in ["Msln", "Upk3b", "Lrrn4", "Dcn", "Krt19", "Cd34", "Pdgfra"]}
print("r3 80 means:", msl, flush=True)
lab = r3.map(ev.S1_best).astype(object)
flag = pd.Series("", index=r3.index, dtype=object)
lab[r3.isin(["29"])] = "cDC2"; flag[r3 == "29"] = "moDC_like"
lab[r3 == "65"] = "ILC1"; lab[r3.isin(["82", "60"])] = "Hepatocyte"
m12 = r3 == "12"; lab[m12] = E.obs.annotation_curated.astype(str)[m12]; flag[m12] = "mixed_r3_12"
lab[r3 == "80"] = "Mesothelial" if (msl["Msln"] > 1 and msl["Upk3b"] > 1) else "Fibroblast"
for c, why in {"1": "hep_LSEC_mixed_nuclei", "61": "hep_LSEC_mixed", "84": "hep_fibro_mixed", "0": "low_UMI_nuclei", "35": "low_UMI_nuclei"}.items():
    lab[r3 == c] = "unresolved"; flag[r3 == c] = why
for c, why in {"39": "KC_with_LSEC_signal", "63": "B_with_myeloid_signal_high_UMI", "67": "KC_with_RBC", "79": "KC_low_UMI",
               "83": "proliferating_single_sample"}.items():
    flag[r3 == c] = why
tdc = sub.index[sub == "3"]; lab[tdc] = "Transitional DC"; flag[tdc] = "tDC_Siglech+Ly6c2+Cd209a+Klk1-"
print(lab.value_counts().to_string(), flush=True)
ORDER = ["Hepatocyte", "Endothelial", "Fibroblast", "Mesothelial", "Cholangiocyte", "KC", "Monocyte", "Neutrophil", "Basophil",
         "cDC1", "cDC2", "Mig. cDCs", "Transitional DC", "pDCs", "B cell", "T cell", "NK cell", "ILC1", "HsPCs", "unresolved"]
assert set(lab) <= set(ORDER), set(lab) - set(ORDER)
cats = [c for c in ORDER if c in set(lab)]
a.obs["annotation_aor"] = pd.Categorical(lab.values, categories=cats); a.obs["aor_flag"] = flag.values
a.obs["leiden_clean_r3"] = r3.values; a.obs["dc_sub"] = sub.reindex(a.obs_names).fillna("").values
a.obsm["X_totalvi"] = E.obsm["X_totalvi"]; a.obsm["X_umap"] = E.obsm["X_umap_totalvi"]
a.obs[["sample", "suspension", "digestion", "sort_gate", "leiden_clean_r3", "dc_sub", "annotation_aor", "aor_flag"]].to_csv(f"{D}/cells_aor.tsv.gz", sep="\t")
pd.crosstab(r3, a.obs.annotation_aor.values).to_csv(f"{D}/r3_labels_aor.tsv", sep="\t")
comp = pd.crosstab(a.obs.annotation_aor, a.obs.suspension); comp["total"] = comp.sum(1); comp["pct"] = (100 * comp.total / comp.total.sum()).round(2)
comp.to_csv(f"{D}/composition_aor.tsv", sep="\t"); print(comp.to_string(), "\ntotal", int(comp.total.sum()), flush=True)
a.uns["note"] = ("ANNOTATION OF RECORD 2026-11-12: strict object (10_strict) + TotalVI 11_totalvi_strict; labels per 13b docstring "
                 "(Table S1 top-50 per r3 cluster + listed overrides). X = clean counts, layers counts = raw UMI.")
for k in a.obs.columns:
    if a.obs[k].dtype == object: a.obs[k] = a.obs[k].astype(str)
a.write_h5ad(f"{D}/guilliams2022_fig1b_aor.h5ad", compression="gzip"); print("written all", a.shape, flush=True)
# figures
pal = dict(zip(cats, sc.pl.palettes.default_20[:len(cats)])); pal["unresolved"] = "#d9d9d9"
n.obs["annotation_aor"] = a.obs.annotation_aor.values; n.obsm["X_umap"] = a.obsm["X_umap"]
fig, axs = plt.subplots(1, 2, figsize=(18, 8))
sc.pl.umap(n, color="annotation_aor", palette=pal, ax=axs[0], show=False, frameon=False, s=1.5, legend_loc="on data", legend_fontsize=7,
           title=f"Guilliams 2022 mouse atlas, rebuilt: {a.n_obs:,} cells + nuclei")
n.obs["suspension"] = a.obs.suspension.values
sc.pl.umap(n, color="suspension", ax=axs[1], show=False, frameon=False, s=1.5, title="suspension")
for ax in axs:
    for c in ax.collections: c.set_rasterized(True)
fig.savefig(f"{D}/figures/umap_aor.png", dpi=200, bbox_inches="tight"); fig.savefig(f"{D}/figures/umap_aor.pdf", dpi=300, bbox_inches="tight"); plt.close(fig)
S1 = "/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/1-s2.0-S0092867421014811-mmc1.xlsx"
x = pd.read_excel(S1, sheet_name=None); dm = {}
for k, v in x.items():
    if k.startswith("Mouse") and k.endswith("DEGs"):
        t = k.replace("Mouse ", "").replace(" DEGs", "")
        if t in cats:
            dm[t] = [g for g in v.sort_values("score", ascending=False).iloc[:, 0].astype(str) if g in n.var_names][:3]
dm["Transitional DC"] = ["Ly6c2", "Cd209a", "Siglech"]; dm["Mesothelial"] = ["Msln", "Upk3b"]
dm = {k: dm[k] for k in cats if k in dm}
k = n[n.obs.annotation_aor != "unresolved"]
k.obs["annotation_aor"] = k.obs.annotation_aor.cat.remove_unused_categories()
dp = sc.pl.dotplot(k, dm, groupby="annotation_aor", standard_scale="var", show=False, return_fig=True)
dp.savefig(f"{D}/figures/dotplot_tableS1_aor.png", dpi=200, bbox_inches="tight"); dp.savefig(f"{D}/figures/dotplot_tableS1_aor.pdf", bbox_inches="tight")
del n, k
import rapids_singlecell as rsc
for susp, tag in (("cell", "cells"), ("nucleus", "nuclei")):
    s = a[a.obs.suspension.astype(str).values == susp].copy()
    t = ad.AnnData(obs=s.obs[[]].copy(), obsm={"Z": s.obsm["X_totalvi"].astype(np.float32)})
    rsc.pp.neighbors(t, use_rep="Z"); rsc.tl.umap(t, random_state=0)
    s.obsm["X_umap_all"] = s.obsm["X_umap"]; s.obsm["X_umap"] = np.asarray(t.obsm["X_umap"])
    if susp == "nucleus":
        for kk in ("adt", "adt_measured"): s.obsm.pop(kk, None)
    s.write_h5ad(f"{D}/guilliams2022_fig1b_aor_{tag}.h5ad", compression="gzip")
    fig, ax = plt.subplots(figsize=(9, 8))
    sc.pl.umap(s, color="annotation_aor", palette=pal, ax=ax, show=False, frameon=False, s=2, legend_loc="on data", legend_fontsize=7,
               title=f"{tag}: {s.n_obs:,}")
    for c in ax.collections: c.set_rasterized(True)
    fig.savefig(f"{D}/figures/umap_aor_{tag}.png", dpi=200, bbox_inches="tight"); plt.close(fig)
    print("written", tag, s.shape, flush=True)
print("done", flush=True)
