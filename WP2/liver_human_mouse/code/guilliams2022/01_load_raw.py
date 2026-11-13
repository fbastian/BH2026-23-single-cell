"""Step 01 - assemble the 46 Fig 1B samples (Guilliams et al. Cell 2022) from GEO GSE192740 into one raw AnnData.

Fig 1A legend: healthy C57Bl/6 mice, 15 ex vivo + 19 in vivo + 12 nuclei samples, 18 CITE-seq (7 ex vivo, 11 in vivo
in the legend; the GEO table gives 8 / 10). GEO rows selected: organism Mouse and title starting with
'Healthy', 'PBS Oral Gavage' or 'Standard diet fed (12 weeks)'.

Input (vault, read-only): GSE192740_RAW.tar (filtered matrices only; 10x h5 for RNA-only samples, dense txt.gz
for CITE-seq samples with gene rows followed by ADT rows), GSE192740_sampleInfo_scRNAseq.tsv.gz,
author repo 4f_TranslateABs.txt (old -> new antibody names).

Output: <OUT>/01_raw/guilliams2022_fig1b_raw.h5ad
  X           raw UMI counts (int32, CSR), var_names = Ensembl IDs, var['symbol']
  obs         sample metadata from GEO + n_counts / n_genes
  obsm['adt'] antibody counts (DataFrame, union of panels, 0 where not measured)
  obsm['adt_measured'] bool DataFrame, True where that antibody was in that sample's panel
Plus <OUT>/01_raw/sample_table.tsv (cells per sample vs GEO 'number of cells').

No filtering, no ambient correction here. Note: the authors ran FastCAR on unfiltered droplets; GEO holds only
filtered matrices, so that step cannot be reproduced exactly.
"""
import gzip
import io
import os
import re
import sys
import tarfile
import tempfile

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

GEO = "/mnt/masld_vault/geo/GSE192740/suppl"
REPO = "/mnt/masld_vault/geo/github-guilliottslab-scripts_GuilliamsEtAll_Cell2022"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/masld_vault/cross_species/guilliams2022"
os.makedirs(f"{OUT}/01_raw", exist_ok=True)

# ---- sample selection -------------------------------------------------------------------------------------------
info = pd.read_csv(f"{GEO}/GSE192740_sampleInfo_scRNAseq.tsv.gz", sep="\t")
info.columns = [c.replace("characteristics: ", "") for c in info.columns]
sel = info[(info["organism"] == "Mouse")
           & info["title"].str.match(r"^(Healthy|PBS Oral Gavage|Standard diet fed \(12 weeks\))")].copy()
counts = sel["digestion method"].value_counts().to_dict()
cite = sel[sel["number of added ABs"] > 0]["digestion method"].value_counts().to_dict()
print("selected", len(sel), counts, "CITE-seq", cite, flush=True)
assert len(sel) == 46 and counts == {"In Vivo": 19, "Ex Vivo": 15, "Nuclei": 12}, "sample selection != Fig 1A"

# ---- antibody name translation ---------------------------------------------------------------------------------
tr = pd.read_csv(f"{REPO}/4f_TranslateABs.txt", sep="\t")
translate = dict(zip(tr["oldName"].astype(str).str.strip(), tr["correctName"].astype(str).str.strip()))

# ---- tar members ------------------------------------------------------------------------------------------------
tar = tarfile.open(f"{GEO}/GSE192740_RAW.tar")
members = {re.sub(r"\.(h5|txt\.gz)$", "", m.name.split("filtered_feature_bc_matrix_")[1]): m for m in tar.getmembers()}

ref_ids = ref_names = None  # gene axis taken from the first h5


def read_h5(m):
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        tmp.write(tar.extractfile(m).read())
        tmp.flush()
        with h5py.File(tmp.name, "r") as h:
            g = h["matrix"]
            ids = [x.decode() for x in g["features/id"][:]]
            names = [x.decode() for x in g["features/name"][:]]
            ftype = [x.decode() for x in g["features/feature_type"][:]]
            bcs = [x.decode() for x in g["barcodes"][:]]
            shape = tuple(g["shape"][:])
            X = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).T.tocsr()
    assert set(ftype) == {"Gene Expression"}, set(ftype)
    return X.astype(np.int32), bcs, ids, names


def read_txt(m):
    df = pd.read_csv(io.BytesIO(tar.extractfile(m).read()), sep="\t", index_col=0, compression="gzip",
                     dtype={0: str}, engine="c")
    df.index = df.index.astype(str)
    return df


sample_rows, blocks, adts = [], [], []
for _, r in sel.iterrows():
    sid = r["shortFileName"]
    m = members[sid]
    gsm = m.name.split("_")[0]
    if m.name.endswith(".h5"):
        X, bcs, ids, names = read_h5(m)
        if ref_ids is None:
            ref_ids, ref_names = ids, names
        assert ids == ref_ids, f"{sid}: gene axis differs from reference h5"
        adt = None
    else:
        if ref_ids is None:
            raise RuntimeError("first selected sample must be h5 to define the gene axis")
        df = read_txt(m)
        g = df.iloc[: len(ref_names)]
        stripped = [re.sub(r"\.\d+$", "", s) if s not in ref_names else s for s in g.index]
        # positional check: names equal after removing R make.unique suffixes
        bad = [(a, b) for a, b in zip(stripped, ref_names) if a != b and re.sub(r"\.\d+$", "", a) != b]
        assert not bad, f"{sid}: gene order mismatch {bad[:5]}"
        X = sp.csr_matrix(g.values.T.astype(np.int32))
        bcs = [b.replace(".", "-") for b in df.columns]
        a = df.iloc[len(ref_names):].T
        a.columns = [translate.get(c, c) for c in a.columns]
        a = a.T.groupby(level=0).sum().T  # merge any names that collapse after translation
        a.index = [f"{sid}_{b}" for b in bcs]
        adt = a.astype(np.int32)
        del df
    obs_names = [f"{sid}_{b}" for b in bcs]
    blocks.append((X, obs_names, sid))
    if adt is not None:
        adts.append(adt)
    n_adt = 0 if adt is None else adt.shape[1]
    sample_rows.append(dict(sample=sid, gsm=gsm, n_cells_file=X.shape[0], n_cells_geo=r["number of cells"],
                            n_abs_geo=r["number of added ABs"], n_adt_rows=n_adt))
    print(f"{sid:8s} {gsm} cells={X.shape[0]:6d} geo={r['number of cells']:6d} adt={n_adt}", flush=True)

# ---- assemble ---------------------------------------------------------------------------------------------------
X = sp.vstack([b[0] for b in blocks], format="csr")
obs_names = [n for b in blocks for n in b[1]]
sample_of = np.concatenate([[b[2]] * b[0].shape[0] for b in blocks])

meta = sel.set_index("shortFileName")
obs = pd.DataFrame(index=pd.Index(obs_names, name="cell"))
obs["sample"] = sample_of
obs["gsm"] = obs["sample"].map(pd.DataFrame(sample_rows).set_index("sample")["gsm"])
obs["title"] = obs["sample"].map(meta["title"])
obs["digestion"] = obs["sample"].map(meta["digestion method"])
obs["chemistry"] = obs["sample"].map(meta["platform"])
obs["n_abs_added"] = obs["sample"].map(meta["number of added ABs"]).astype(int)
obs["suspension"] = np.where(obs["digestion"] == "Nuclei", "nucleus", "cell")
obs["mouse"] = obs["title"].str.extract(r"_(Mouse_\d+)_")[0]
obs["condition"] = obs["title"].str.replace(r"_Mouse_\d+_Sample_\d+$", "", regex=True)
obs["cite_seq"] = obs["n_abs_added"] > 0
obs["n_counts"] = np.asarray(X.sum(1)).ravel()
obs["n_genes"] = np.diff(X.indptr)
for c in ["sample", "gsm", "title", "digestion", "chemistry", "suspension", "mouse", "condition"]:
    obs[c] = obs[c].astype("category")

var = pd.DataFrame(index=pd.Index(ref_ids, name="gene_id"))
var["symbol"] = ref_names
adata = ad.AnnData(X=X, obs=obs, var=var)

adt_all = pd.concat(adts, axis=0, join="outer")
measured = adt_all.notna()
adt_all = adt_all.reindex(adata.obs_names).fillna(0).astype(np.int32)
measured = measured.reindex(adata.obs_names).fillna(False).astype(bool)
adata.obsm["adt"] = adt_all
adata.obsm["adt_measured"] = measured
adata.uns["provenance"] = dict(
    source="GEO GSE192740 filtered matrices; Guilliams et al. Cell 2022 Fig 1B sample set",
    selection="organism Mouse & title ~ ^(Healthy|PBS Oral Gavage|Standard diet fed (12 weeks))",
    adt_names="translated with author repo 4f_TranslateABs.txt",
    note="no QC, no ambient correction (FastCAR needs unfiltered droplets, not on GEO)",
)

st = pd.DataFrame(sample_rows)
st["diff_file_minus_geo"] = st["n_cells_file"] - st["n_cells_geo"]
st.to_csv(f"{OUT}/01_raw/sample_table.tsv", sep="\t", index=False)
print(st.to_string(index=False))
print("TOTAL", adata.shape, "ADT panel union", adt_all.shape[1], flush=True)
adata.write_h5ad(f"{OUT}/01_raw/guilliams2022_fig1b_raw.h5ad", compression="gzip")
print("written", flush=True)
