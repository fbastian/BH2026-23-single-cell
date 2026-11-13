"""Step 13a - annotation of record, evidence tables (user go 2026-11-12 16:30 CET).
Input: 10_strict (234,552, X = clean counts) + 11_totalvi_strict embeddings (leiden_clean_r3, 85 clusters).
1. Table S1 (mmc1) top-50 mouse DEGs per Fig 1B type, cluster-level on r3, z within cell- / nucleus-majority clusters (as 06f).
2. Classic marker means per r3 cluster, curated-label mode, doublet score, UMIs, suspension, digestion, sort gate.
3. DC block (r3 clusters whose curated mode is pDCs or cDC2): Leiden r1 on the TotalVI latent of that subset, marker profile
   per subcluster, for the transitional DC label (Siglech+ Ly6c2+ Cd209a+, Klk1-; check in 12b).
Writes <OUT>/13_annotation_of_record/: r3_evidence.tsv, r3_tableS1_scores.tsv, dc_sub_evidence.tsv, dc_sub.tsv.gz,
lognorm cache not written (recomputed in 13b)."""
import os
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
OUT = "/mnt/masld_vault/cross_species/guilliams2022"; D = f"{OUT}/13_annotation_of_record"; os.makedirs(f"{D}/figures", exist_ok=True)
S1 = "/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/1-s2.0-S0092867421014811-mmc1.xlsx"
E = ad.read_h5ad(f"{OUT}/11_totalvi_strict/embeddings.h5ad")
a = ad.read_h5ad(f"{OUT}/10_strict/guilliams2022_fig1b_strict.h5ad")
assert (a.obs_names == E.obs_names).all()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sym = a.var.symbol.astype(str).values
cl = E.obs.leiden_clean_r3.astype(str).values; cats = sorted(set(cl), key=int)
codes = pd.Categorical(cl, categories=cats).codes
ind = sp.csr_matrix((np.ones(len(cl)), (codes, np.arange(len(cl)))), shape=(len(cats), len(cl))); ncl = np.asarray(ind.sum(1)).ravel()
def means(genes, X=a.X, ind=ind, ncl=ncl, idx=cats):
    cols = np.where(np.isin(sym, genes))[0]
    return pd.DataFrame(np.asarray((ind @ X[:, cols]).todense()) / ncl[:, None], index=idx, columns=sym[cols]).T.groupby(level=0).mean().T
x = pd.read_excel(S1, sheet_name=None); markers = {}
for k, v in x.items():
    if k.startswith("Mouse") and k.endswith("DEGs"):
        g = v.sort_values("score", ascending=False).iloc[:50, 0].astype(str)
        markers[k.replace("Mouse ", "").replace(" DEGs", "")] = [s for s in g if s in set(sym)]
G = means(sorted({s for v in markers.values() for s in v}))
o = E.obs
fnuc = pd.Series(o.suspension.astype(str).values == "nucleus").groupby(cl).mean().reindex(cats)
grp = np.where(fnuc.values > 0.5, "nucleus", "cell")
score = pd.DataFrame(index=cats, columns=list(markers), dtype=float)
for gname in ("cell", "nucleus"):
    sel = grp == gname; Gs = G.loc[sel]; Z = (Gs - Gs.mean()) / Gs.std().replace(0, np.nan)
    for t_, gl in markers.items():
        score.loc[Gs.index, t_] = Z[gl].mean(1)
score.to_csv(f"{D}/r3_tableS1_scores.tsv", sep="\t", float_format="%.3g")
srt = np.sort(score.values, 1)
CLASSIC = ["Alb", "Cyp2e1", "Cyp2f2", "Clec4f", "Vsig4", "Stab2", "Kdr", "Ptprc", "Cd3e", "Cd8a", "Nkg7", "Ncr1", "Gzma", "Cd79a",
           "Siglech", "Klk1", "Xcr1", "Cd209a", "Ccr7", "S100a8", "Ly6c2", "Lyz2", "Reln", "Dcn", "Krt19", "Mcpt8", "Kit", "Cd34",
           "Mki67", "Hba-a1"]
C = means(CLASSIC)[[g for g in CLASSIC if g in sym]]
gb = o.groupby(cl, observed=True)
ev = pd.DataFrame({"n": ncl.astype(int), "frac_nuc": fnuc.values, "S1_best": score.idxmax(1), "S1_score": srt[:, -1],
                   "S1_second": score.apply(lambda r: r.drop(r.idxmax()).idxmax(), axis=1), "margin": srt[:, -1] - srt[:, -2]}, index=cats)
ev["curated_mode"] = gb.annotation_curated.agg(lambda s: s.astype(str).mode().iat[0]).reindex(cats)
ev["curated_frac"] = gb.annotation_curated.agg(lambda s: (s.astype(str) == s.astype(str).mode().iat[0]).mean()).reindex(cats)
ev["flag"] = gb.annotation_flag.agg(lambda s: ",".join(sorted(set(s.astype(str)) - {"", "nan", "none"}))).reindex(cats)
ev["dbl"] = gb.scDblFinder_score.median().reindex(cats); ev["umis"] = gb.umis_final.median().reindex(cats)
ev["exvivo"] = gb.digestion.agg(lambda s: (s.astype(str) == "Ex Vivo").mean()).reindex(cats)
ev["n_samp"] = gb["sample"].nunique().reindex(cats)
ev["top_samp"] = gb["sample"].agg(lambda s: f"{s.astype(str).mode().iat[0]}:{(s.astype(str) == s.astype(str).mode().iat[0]).mean():.2f}").reindex(cats)
ev = ev.join(C.round(2))
ev.to_csv(f"{D}/r3_evidence.tsv", sep="\t", float_format="%.3g")
pd.set_option("display.width", 500); pd.set_option("display.max_columns", 80)
print(ev.drop(columns=CLASSIC, errors="ignore").sort_values(["S1_best", "n"]).to_string(float_format=lambda v: f"{v:.2f}"), flush=True)
print(C.round(1).loc[ev.sort_values(["S1_best", "n"]).index].assign(lab=ev.S1_best).to_string(), flush=True)
# DC block subclustering
import rapids_singlecell as rsc
dcl = ev.index[ev.curated_mode.isin(["pDCs", "cDC2"])]
m = np.isin(cl, dcl)
t = ad.AnnData(obs=o[m][["sample", "digestion", "annotation_curated"]].copy(), obsm={"Z": E.obsm["X_totalvi"][m].astype(np.float32)})
rsc.pp.neighbors(t, use_rep="Z"); rsc.tl.leiden(t, resolution=1.0, random_state=0, key_added="sub")
sub = t.obs["sub"].astype(str).values; scats = sorted(set(sub), key=int)
sind = sp.csr_matrix((np.ones(m.sum()), (pd.Categorical(sub, categories=scats).codes, np.arange(m.sum()))), shape=(len(scats), m.sum()))
sn = np.asarray(sind.sum(1)).ravel()
DCG = ["Siglech", "Ccr9", "Klk1", "Ly6d", "Tcf4", "Ly6c2", "Zeb2", "Cx3cr1", "Cd209a", "Itgax", "Sirpa", "Clec4a4", "Mgl2", "H2-Ab1",
       "Esam", "Lyz2", "Mki67", "Ccr2"]
X = a.X[np.flatnonzero(m)]
cols = np.where(np.isin(sym, DCG))[0]
SM = pd.DataFrame(np.asarray((sind @ X[:, cols]).todense()) / sn[:, None], index=scats, columns=sym[cols]).T.groupby(level=0).mean().T[[g for g in DCG if g in sym]]
og = t.obs.groupby(sub)
sev = pd.DataFrame({"n": sn.astype(int), "r3_mode": pd.Series(cl[m]).groupby(sub).agg(lambda s: s.mode().iat[0]).reindex(scats).values,
                    "curated_mode": og.annotation_curated.agg(lambda s: s.astype(str).mode().iat[0]).reindex(scats).values,
                    "exvivo": og.digestion.agg(lambda s: (s.astype(str) == "Ex Vivo").mean()).reindex(scats).values,
                    "n_samp": og["sample"].nunique().reindex(scats).values,
                    "umis": pd.Series(o.umis_final.values[m]).groupby(sub).median().reindex(scats).values,
                    "dbl": pd.Series(o.scDblFinder_score.values[m]).groupby(sub).median().reindex(scats).values}, index=scats).join(SM.round(2))
U = E.obsm["X_umap_totalvi"][m]
sev["umap"] = [f"{U[sub == s, 0].mean():.1f},{U[sub == s, 1].mean():.1f}" for s in scats]
sev.to_csv(f"{D}/dc_sub_evidence.tsv", sep="\t", float_format="%.3g")
print(sev.to_string(float_format=lambda v: f"{v:.2f}"), flush=True)
pd.DataFrame({"dc_sub": sub}, index=t.obs_names).to_csv(f"{D}/dc_sub.tsv.gz", sep="\t")
print("done", flush=True)
