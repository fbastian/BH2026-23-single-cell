"""HLiCA figure data (2026-11-12): per-cell UMAP (deposited X_umap) + labels for all 524,699 cells; lineage dot plot statistics
(HLiCA canonical lineage markers, Zenodo 20380927) on the SAMap subsample (50,400 cells, <= 1,200 per author_cell_type, raw counts ->
normalise_total 1e4 + log1p), grouped by HLiCA lineage (colour_palette.R alpha). Outputs <D>/data/hlica_*."""
import sys; sys.path.insert(0, "/home/ubuntu/liver_cross_species/figures/code")
import anndata as ad, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from common import D, HLICA_LIN
h = ad.read_h5ad("/mnt/masld_vault/geo/cxg-059202e1-1f1b-483f-9151-f3a25a380c39/HLiCA_human_liver_all_cells.h5ad", backed="r")
o = h.obs
pd.DataFrame({"label": o.author_cell_type.astype(str).values, "suspension": o.suspension_type.astype(str).values, "study": o.STUDY.astype(str).values,
              "donor": o.donor_id.astype(str).values, "u1": h.obsm["X_umap"][:, 0], "u2": h.obsm["X_umap"][:, 1]}, index=o.index).to_parquet(f"{D}/data/hlica_umap.parquet")
mk = pd.read_csv("/mnt/masld_vault/geo/zenodo-20380927/HLiCA_lineage_markers.csv")
print(mk.groupby(["Lineage", "Marker Type"]).Gene.apply(list).to_string())
s = ad.read_h5ad("/mnt/masld_vault/cross_species/samap_run/input/hs_samap_ready.h5ad")
sc.pp.normalize_total(s, target_sum=1e4); sc.pp.log1p(s)
genes = [g for g in mk.Gene if g in s.var_names]
lin = s.obs.author_cell_type.astype(str).map(HLICA_LIN).values
X = s[:, genes].X; X = X.toarray() if sp.issparse(X) else X
pd.DataFrame(X, columns=genes).groupby(lin).mean().to_csv(f"{D}/data/hlica_dot_mean.tsv", sep="\t")
pd.DataFrame(X > 0, columns=genes).groupby(lin).mean().to_csv(f"{D}/data/hlica_dot_frac.tsv", sep="\t")
mk[mk.Gene.isin(genes)].to_csv(f"{D}/data/hlica_dot_genes.tsv", sep="\t", index=False)
print("done")
