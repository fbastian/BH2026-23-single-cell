"""Step 6 - sanity check of the MMseqs2 maps against MGI curated 1:1 orthologs (HOM_MouseHumanSequence.rpt, fetched 2026-11-12;
Ensembl BioMart and HCOP downloads failed from the VM). 1:1 = homology class with exactly one mouse and one human gene.
Restricted to human protein-coding genes (HLiCA var feature_type). Per direction: of the 1:1 pairs whose two genes are in our
gene tables (matched by symbol = samap_gene), the fraction where the ortholog is (a) the top hit by bits, (b) among the hits at all.
Output refs/orthologs/ortholog_check.tsv + printed summary."""
import anndata as ad, pandas as pd
R = "/mnt/masld_vault/cross_species/samap_maps/liver_hsmm"
cols = ["q", "t", "pident", "alnlen", "mm", "gap", "qs", "qe", "ts", "te", "evalue", "bits"]
mgi = pd.read_csv(f"{R}/refs/orthologs/HOM_MouseHumanSequence.rpt", sep="\t", dtype=str)
g = mgi.groupby("DB Class Key")
one = g.filter(lambda x: (x["NCBI Taxon ID"] == "10090").sum() == 1 and (x["NCBI Taxon ID"] == "9606").sum() == 1)
p = one.pivot_table(index="DB Class Key", columns="NCBI Taxon ID", values="Symbol", aggfunc="first").rename(columns={"9606": "hs", "10090": "mm"})
hv = ad.read_h5ad("/mnt/masld_vault/geo/cxg-059202e1-1f1b-483f-9151-f3a25a380c39/HLiCA_human_liver_all_cells.h5ad", backed="r").var
pc = set(hv.feature_name[hv.feature_type.astype(str) == "protein_coding"].astype(str))
hsn = set(pd.read_csv(f"{R}/refs/mappings/hs.gene_names.tsv", sep="\t").samap_gene); mmn = set(pd.read_csv(f"{R}/refs/mappings/mm.gene_names.tsv", sep="\t").samap_gene)
p = p[p.hs.isin(pc) & p.hs.isin(hsn) & p.mm.isin(mmn)]
print("MGI 1:1 pairs, human protein-coding, both genes in our tables:", len(p), flush=True)
out = []
for d, qcol, tcol in (("hs_to_mm", "hs", "mm"), ("mm_to_hs", "mm", "hs")):
    m = pd.read_csv(f"{R}/maps_mmseqs/hsmm/{d}.txt", sep="\t", header=None, names=cols, usecols=["q", "t", "bits", "pident"])
    m = m[m.q.isin(set(p[qcol]))]
    top = m.sort_values("bits", ascending=False).drop_duplicates("q").set_index("q").t
    pairs = set(zip(m.q, m.t))
    exp = p.set_index(qcol)[tcol]
    hit = exp.index.isin(top.index)
    r = {"direction": d, "n_pairs": len(exp), "query_has_hits": hit.mean(),
         "ortholog_top_hit": (top.reindex(exp.index) == exp).mean(),
         "ortholog_among_hits": pd.Series([(q, t) in pairs for q, t in exp.items()]).mean()}
    out.append(r); print(r, flush=True)
pd.DataFrame(out).to_csv(f"{R}/refs/orthologs/ortholog_check.tsv", sep="\t", index=False, float_format="%.4f")
