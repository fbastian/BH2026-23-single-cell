"""Step 12b - profile the cDC2 -> pDC streak along its axis (bins of the projection between the cDC2 and pDC UMAP centres).
Doublet test: doublets of pDC x cDC2 would carry ~2x UMIs, higher scDblFinder score and co-expression of both full
programmes; transitional DCs (tDC) carry a graded Siglech/Ly6c2 -> Cd209a/Itgax switch with Zeb2 and normal UMIs.
Output <OUT>/11_totalvi_strict/dc_streak/streak_profile.tsv, streak_cells.tsv.gz, streak_samples.tsv."""
import anndata as ad, numpy as np, pandas as pd, scanpy as sc
OUT = "/mnt/masld_vault/cross_species/guilliams2022"; D = f"{OUT}/11_totalvi_strict/dc_streak"
E = ad.read_h5ad(f"{OUT}/11_totalvi_strict/embeddings.h5ad")
U = E.obsm["X_umap_totalvi"]; lab = E.obs.annotation_curated.astype(str).values
# centres = medians of the core clusters (pDC r3 30, cDC2 r3 2)
r3 = E.obs.leiden_clean_r3.astype(str).values
c0, c1 = np.median(U[r3 == "2"], 0), np.median(U[r3 == "30"], 0)
v = c1 - c0; L = np.linalg.norm(v); v = v / L
t = (U - c0) @ v / L; d = np.abs((U - c0) @ np.array([-v[1], v[0]]))
sel = (t > -0.1) & (t < 1.1) & (d < 0.6)
A = ad.read_h5ad(f"{OUT}/10_strict/guilliams2022_fig1b_strict.h5ad", backed="r")
m = A[np.flatnonzero(sel)].to_memory(); s = ad.AnnData(X=m.X.copy(), obs=E.obs.iloc[np.flatnonzero(sel)].copy(), var=A.var.copy()); del m
sc.pp.normalize_total(s, target_sum=1e4); sc.pp.log1p(s)
s.var_names = s.var.symbol.astype(str).values; s.var_names_make_unique()
g = ["Siglech", "Ly6c2", "Ccr9", "Klk1", "Tcf4", "Cd209a", "Itgax", "Sirpa", "Zeb2", "Cx3cr1", "Lyz2", "H2-Ab1", "Esam", "Clec4a4", "Mki67"]
X = pd.DataFrame(s[:, g].X.toarray(), index=s.obs_names, columns=g)
o = s.obs; o["t"] = t[sel]; o["bin"] = pd.cut(o.t, [-0.1, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.1])
o["co_pDC_cDC2"] = (X.Siglech > 1) & (X.Cd209a > 1)
grp = o.groupby("bin", observed=True)
P = pd.DataFrame({"n": grp.size(), "pDC_lab": grp.annotation_curated.agg(lambda x: (x.astype(str) == "pDCs").mean()),
                  "cDC2_lab": grp.annotation_curated.agg(lambda x: (x.astype(str) == "cDC2").mean()),
                  "umis_med": grp.umis_final.median(), "dbl_med": grp.scDblFinder_score.median(),
                  "dbl_called_raw": grp.scDblFinder_class.agg(lambda x: (x.astype(str) == "doublet").mean()),
                  "co_Siglech_Cd209a": grp.co_pDC_cDC2.mean(), "exvivo": grp.digestion.agg(lambda x: (x.astype(str) == "Ex Vivo").mean()),
                  "n_samples": grp["sample"].nunique()})
P = P.join(X.groupby(o.bin.values, observed=True).mean())
pd.set_option("display.width", 400)
print(P.to_string(float_format=lambda v: f"{v:.2f}"))
P.to_csv(f"{D}/streak_profile.tsv", sep="\t", float_format="%.3g")
mid = o.t.between(0.25, 0.85)
st = pd.crosstab(o.loc[mid, "sample"].astype(str), o.loc[mid, "sort_gate"].astype(str))
st["n"] = st.sum(1); st = st.sort_values("n", ascending=False)
print("\nstreak middle (t 0.25-0.85):", int(mid.sum()), "cells, samples:\n", st.head(15).to_string())
st.to_csv(f"{D}/streak_samples.tsv", sep="\t")
pd.concat([o[["sample", "sort_gate", "digestion", "annotation_curated", "leiden_clean_r3", "t", "bin", "umis_final", "scDblFinder_score"]].astype({"bin": str}), X], axis=1).to_csv(f"{D}/streak_cells.tsv.gz", sep="\t")
