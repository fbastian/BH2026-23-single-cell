"""SAMap run HLiCA (hs) x Guilliams 2022 mouse (mm), sc-samap 3.0.1 (Project 23 WP2 recipe: SAMAP(filenames, f_maps), sm.run(),
get_mapping_scores n_top=0, GenePairFinder align_thr 0.1). Maps: our liver MMseqs2 maps (samap_maps/liver_hsmm/maps_mmseqs/hsmm).
Keys: hs author_cell_type (48 HLiCA labels), mm annotation_aor. Outputs /mnt/masld_vault/cross_species/samap_run/out/."""
import os, pickle, time
import numpy as np, pandas as pd
from scipy import sparse
from samap import SAMAP
from samap.analysis import get_mapping_scores, GenePairFinder
from samap.utils import save_samap
t0 = time.time(); log = lambda s: print(f"[{(time.time() - t0) / 60:6.1f} min] {s}", flush=True)
I = "/mnt/masld_vault/cross_species/samap_run/input"; O = "/mnt/masld_vault/cross_species/samap_run/out"
MAPS = "/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/maps_mmseqs/"
sm = SAMAP({"hs": f"{I}/hs_samap_ready.h5ad", "mm": f"{I}/mm_samap_ready.h5ad"}, f_maps=MAPS)
log("SAMAP object built")
sm.run(ncpus=16)
log("sm.run done")
keys = {"hs": "author_cell_type", "mm": "annotation_aor"}
D, MT = get_mapping_scores(sm, keys, n_top=0)
pd.DataFrame(D).to_csv(f"{O}/mapping_D.tsv", sep="\t"); MT.to_csv(f"{O}/mapping_table.tsv", sep="\t")
log("mapping scores written")
sm.samap.adata.write_h5ad(f"{O}/samap_joint.h5ad")
for sp_, sam in sm.sams.items():
    X = sam.adata.X
    if sparse.issparse(X):
        X = X.tocsr(copy=False); X.indices = X.indices.astype(np.int64, copy=False); X.indptr = X.indptr.astype(np.int64, copy=False); sam.adata.X = X
try:
    gpf = GenePairFinder(sm, keys=keys)
    gp = gpf.find_all(align_thr=0.1)
    pd.DataFrame(gp).to_csv(f"{O}/gene_pairs_thr0.1.tsv", sep="\t", index=False)
    log("gene pairs written")
except Exception as e:
    log(f"GenePairFinder failed: {e!r}")
save_samap(sm, f"{O}/samap_object.pkl")
log("done")
