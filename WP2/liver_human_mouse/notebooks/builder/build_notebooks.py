"""Build the reproduction notebooks (2026-11-12). Run: ~/miniforge3/envs/masld_atlas/bin/python build_notebooks.py
Each step's script is read verbatim from the project code folders; only the root /mnt/masld_vault/cross_species/ is
replaced by the REPRO mirror root and the figures code folder by WORK. Light notebooks first (01-08, minutes each, they read
the outputs of record of heavy steps), heavy notebooks last (09-12, hours, guarded by RUN_HEAVY)."""
import os
import nbformat as nbf

H = os.path.expanduser("~/liver_cross_species")
NB = f"{H}/notebooks"; WORK = f"{NB}/work"
REC, REP = "/mnt/masld_vault/cross_species/", "/mnt/masld_vault/cross_species_repro/"
G, SM, SR, FG = f"{H}/guilliams2022/code", f"{H}/samap_maps/code", f"{H}/samap_run/code", f"{H}/figures/code"
KERNEL = {"name": "masld_atlas", "display_name": "masld_atlas", "language": "python"}


def src(path):
    t = open(path).read().replace(REC, REP).replace(f"{H}/figures/code", WORK)
    return t.replace("~/liver_cross_species/samap_maps", f"{NB}/work/samap_maps_run")


def md(t): return nbf.v4.new_markdown_cell(t.strip())
def code(t): return nbf.v4.new_code_cell(t.strip())


def step(path, run_line, title=None, note=None):
    name = os.path.basename(path)
    cells = []
    if title: cells.append(md(f"### {title}" + (f"\n\n{note.strip()}" if note else "")))
    cells.append(code(f"%%writefile {WORK}/{name}\n" + src(path)))
    if run_line: cells.append(code(run_line))
    return cells


SETUP = code(f"""%run -i {NB}/repro_setup.py
print("RECORD", RECORD); print("REPRO ", REPRO); print("WORK  ", WORK)""")


def header(title, body, runtime, envs):
    return md(f"""# {title}

{body.strip()}

**Runtime:** {runtime}
**Envs used (each step runs as a subprocess in its own env):** {envs}
**Kernel:** `masld_atlas` (orchestration only).

Outputs are written under `/mnt/masld_vault/cross_species_repro/` (mirror of `/mnt/masld_vault/cross_species/`, the outputs
of record). Inputs this notebook does not rebuild are linked from the record; see `repro_setup.py` for the safety rules.""")


def save(name, cells):
    nb = nbf.v4.new_notebook(); nb.cells = cells
    nb.metadata["kernelspec"] = KERNEL; nb.metadata["language_info"] = {"name": "python"}
    nbf.write(nb, f"{NB}/{name}"); print("written", name, len(cells), "cells")


# ---------------------------------------------------------------- 01
save("01_raw_object.ipynb", [
    header("01 - Raw object: the 46 Fig 1B samples from GEO",
           """Guilliams et al., Cell 2022 (doi 10.1016/j.cell.2021.12.018), mouse healthy atlas, Fig 1B (paper: 185,894 cells and
nuclei after QC). GEO GSE192740 holds per-sample filtered matrices only: 46 samples (15 ex vivo, 19 in vivo, 12 nuclei; 18 with
CITE-seq). Result of record: 389,056 cells, 31,053 genes, 173 ADTs.""", "~10-15 min (reads the GEO tar)", "masld_integration"),
    SETUP,
    code('guard(30)\nfresh("guilliams2022/01_raw")'),
    *step(f"{G}/01_load_raw.py", 'run("01_load_raw.py", "masld_integration")'),
    md("### Check against the record"),
    code('compare_tsv("guilliams2022/01_raw/sample_table.tsv", index_col=None)'),
])

# ---------------------------------------------------------------- 02
save("02_qc_authors_rule.ipynb", [
    header("02 - QC: the authors' rule + Scrublet",
           """Authors' 3a_SCT.R: per sample, 3 MADs on log genes and log UMI (both sides) and % mito (upper). Their hand removal of doublet
clusters is stood in for by Scrublet here (per sample, GPU). Result of record: 328,041 of 389,056 kept.""", "~5 min (GPU)", "rapids-sc"),
    SETUP,
    code('guard(20)\nlink("guilliams2022/01_raw")\nfresh("guilliams2022/02_qc")'),
    *step(f"{G}/02_qc.py", 'run("02_qc.py", "rapids-sc", gpu=0)'),
    md("### Check against the record"),
    code('compare_tsv("guilliams2022/02_qc/qc_per_sample.tsv")'),
])

# ---------------------------------------------------------------- 03
save("03_draft_annotation.ipynb", [
    header("03 - Draft annotation from Table S1",
           """Leiden r2/r3 on the draft TotalVI latent (heavy step, notebook 12; the record's latent is linked here), then cluster-level
scoring with the paper's Table S1 (mmc1) top-50 mouse DEGs for the 17 Fig 1B types. Leiden on GPU: if the clusters differ from the
record, the record's clusters are restored, because later hand decisions (doublet clusters in 04, overrides in 05) name record
cluster numbers.""", "~5 min (GPU)", "rapids-sc"),
    SETUP,
    code('link("guilliams2022/02_qc", "guilliams2022/03_integration")\nfresh("guilliams2022/05_annotation_draft")'),
    *step(f"{G}/05_draft_annotation.py", 'run("05_draft_annotation.py", "rapids-sc", gpu=0)'),
    md("### Check against the record (and restore record clusters if the GPU Leiden differs)"),
    code('''import pandas as pd, shutil
a = pd.read_parquet(rec("guilliams2022/05_annotation_draft/cells_annotation_draft.parquet"))
b = pd.read_parquet(rep("guilliams2022/05_annotation_draft/cells_annotation_draft.parquet"))
same = {c: float((a[c].astype(str) == b[c].astype(str)).mean()) for c in ["leiden_totalvi_r2", "leiden_totalvi_r3", "annotation_draft"]}
print(same)
if same["leiden_totalvi_r3"] < 1:
    for f in ["cells_annotation_draft.parquet", "cluster_labels.tsv", "cluster_scores.tsv", "composition_draft.tsv"]:
        shutil.copy2(rec(f"guilliams2022/05_annotation_draft/{f}"), rep(f"guilliams2022/05_annotation_draft/{f}"))
    print("GPU Leiden differs from the record: record clusters and labels restored (downstream hand decisions use them)")'''),
    code('show("guilliams2022/05_annotation_draft/figures/umap_annotation_draft.png")'),
])

# ---------------------------------------------------------------- 04
save("04_doublets_ambient_clean.ipynb", [
    header("04 - Doublets, ambient RNA and the clean object",
           """Assembles the per-sample scDblFinder + decontX results (heavy, notebook 10; record linked) and builds the clean object.

Decisions (2026-11-11):
* decontX strength from the ladder test (06d): fixed prior delta (1, 300) for the 13 in vivo cell samples with non-epithelial
  Alb > 5 CP10k; raw counts for ex vivo, clean in vivo samples and all nuclei; CS140 raw + flag `ambient_dominated`.
* Scrublet calls were already removed by the draft QC, so the MASLD consensus rule (Scrublet AND scDblFinder) cannot apply;
  the 8 TotalVI r3 clusters with > 65 % scDblFinder doublets are removed (6,982 cells), as the authors removed doublet clusters.
* Re-annotation with Table S1, z-scores within cell-majority and nucleus-majority clusters separately.

Result of record: 321,059 cells and nuclei (`07_clean`).""", "~25 min", "rapids-sc"),
    SETUP,
    code('''guard(30)
link("guilliams2022/02_qc", "guilliams2022/03_integration", "guilliams2022/05_annotation_draft")
fresh("guilliams2022/06_doublets_ambient")
link("guilliams2022/06_doublets_ambient/per_sample", "guilliams2022/06_doublets_ambient/ladder", "guilliams2022/06_doublets_ambient/final_dx")'''),
    md("### decontX ladder (record of notebook 10): Alb removed from non-epithelial cells vs Alb kept in hepatocytes"),
    code('''import pandas as pd, glob
pd.concat([pd.read_csv(f, sep="\\t").assign(sample=f.split("/")[-1].split(".")[0]) for f in sorted(glob.glob(rep("guilliams2022/06_doublets_ambient/ladder/*.ladder.tsv")))])'''),
    *step(f"{G}/06c_assemble.py", 'run("06c_assemble.py", "rapids-sc")', "06c - assemble per-sample results"),
    code('compare_tsv("guilliams2022/06_doublets_ambient/per_sample_summary.tsv")'),
    code('fresh("guilliams2022/07_clean")'),
    *step(f"{G}/06f_clean_annotate_split.py", 'run("06f_clean_annotate_split.py", "rapids-sc", gpu=0)', "06f - clean object, re-annotation, cells / nuclei split"),
    code('compare_tsv("guilliams2022/07_clean/cluster_labels_clean.tsv")\ncompare_tsv("guilliams2022/07_clean/composition_clean.tsv")'),
])

# ---------------------------------------------------------------- 05
save("05_strict_qc.ipynb", [
    header("05 - Manual overrides, streak diagnosis and strict QC",
           """* 07: hand overrides of 4 nuclei clusters (`07_clean/nuclei_overrides.tsv`: 4 -> Hepatocyte flagged possible doublet,
  76 -> Cholangiocyte, 0 and 62 -> unresolved).
* 09: why the UMAP still had streaks after the clean-up (needs the clean TotalVI latent, heavy notebook 12; record linked).
* 10: strict QC approximating the authors' "strict QC controls eliminating these potential doublets" (Cell 2022 p.382):
  per-cell scDblFinder doublets, < 200 UMI after decontX, between-type (bridge) cells, FACS sort-gate mismatches.

Result of record: 234,552 kept (sort mismatch 39,843; scDblFinder 35,042; < 200 UMI 8,522; bridge 3,100).""", "~20 min", "rapids-sc"),
    SETUP,
    md("### 07 - overrides (rewrites the 07_clean files in place, so it needs notebook 04's rebuilt 07_clean)"),
    code('''import shutil
if is_real("guilliams2022/07_clean"):
    shutil.copy2(rec("guilliams2022/07_clean/nuclei_overrides.tsv"), rep("guilliams2022/07_clean/nuclei_overrides.tsv"))
    print(open(rep("guilliams2022/07_clean/nuclei_overrides.tsv")).read())
    DO_07 = True
else:
    link("guilliams2022/07_clean"); DO_07 = False
    print("07_clean not rebuilt (run notebook 04 first): the record's 07_clean already carries the overrides; step 07 skipped")'''),
    *step(f"{G}/07_apply_overrides.py", 'if DO_07: run("07_apply_overrides.py", "rapids-sc")'),
    md("### 09 - streak diagnostic"),
    code('link("guilliams2022/08_totalvi_clean")   # clean TotalVI latent (heavy notebook 12) from the record unless rebuilt\nfresh("guilliams2022/08_totalvi_clean/streaks")'),
    *step(f"{G}/09_streak_diagnostic.py", 'run("09_streak_diagnostic.py", "rapids-sc", gpu=0)'),
    code('show("guilliams2022/08_totalvi_clean/streaks/umap_streak_drivers.png")'),
    md("### 10 - strict QC"),
    code('guard(30)\nfresh("guilliams2022/10_strict")'),
    *step(f"{G}/10_strict_qc.py", 'run("10_strict_qc.py", "rapids-sc", gpu=0)'),
    code('compare_tsv("guilliams2022/10_strict/per_sample_strict.tsv")\ncompare_tsv("guilliams2022/10_strict/removal_by_type.tsv")'),
])

# ---------------------------------------------------------------- 06
save("06_annotation_of_record.ipynb", [
    header("06 - The cDC2-pDC streak and the annotation of record",
           """* 12 / 12b: the streak between cDC2 and pDCs on the strict TotalVI UMAP (latent from heavy notebook 12; record linked).
  Not doublets (UMIs lower, not higher; scDblFinder score flat; Klk1 absent; Ly6c2 peaks in the middle): transitional DCs.
* 13a: evidence per strict Leiden r3 cluster (Table S1 scores, classic markers) + DC-block subclusters.
* 13b: the annotation of record (overrides listed in the 13b docstring). Result of record: 16 of 17 Fig 1B types (no HsPC
  cluster) + Transitional DC (735) + Mesothelial (359); 6,859 unresolved (kept, excluded from SAMap).

The DC subclusters come from a GPU Leiden; 13b names subcluster 3, so the record's subclusters are restored if they differ.""",
           "~25 min", "masld_integration (12, 12b), rapids-sc (13a, 13b)"),
    SETUP,
    code('link("guilliams2022/10_strict", "guilliams2022/11_totalvi_strict")   # strict TotalVI latent from the record unless rebuilt\nfresh("guilliams2022/11_totalvi_strict/dc_streak")'),
    *step(f"{G}/12_dc_streak_check.py", 'run("12_dc_streak_check.py", "masld_integration")', "12 - zoom on the DC region"),
    *step(f"{G}/12b_dc_streak_profile.py", 'run("12b_dc_streak_profile.py", "masld_integration")', "12b - profile along the streak"),
    code('compare_tsv("guilliams2022/11_totalvi_strict/dc_streak/streak_profile.tsv")\nshow("guilliams2022/11_totalvi_strict/dc_streak/umap_dc_zoom.png")'),
    code('guard(20)\nfresh("guilliams2022/13_annotation_of_record")'),
    *step(f"{G}/13a_aor_scores.py", 'run("13a_aor_scores.py", "rapids-sc", gpu=0, tail=10)', "13a - evidence tables"),
    code('''import pandas as pd, shutil
a = pd.read_csv(rec("guilliams2022/13_annotation_of_record/dc_sub.tsv.gz"), sep="\\t", index_col=0)
b = pd.read_csv(rep("guilliams2022/13_annotation_of_record/dc_sub.tsv.gz"), sep="\\t", index_col=0)
same = float((a.dc_sub.astype(str) == b.dc_sub.reindex(a.index).astype(str)).mean()); print("DC subclusters identical to record:", same)
if same < 1:
    for f in ["dc_sub.tsv.gz", "dc_sub_evidence.tsv"]:
        shutil.copy2(rec(f"guilliams2022/13_annotation_of_record/{f}"), rep(f"guilliams2022/13_annotation_of_record/{f}"))
    print("record DC subclusters restored (13b names subcluster 3)")
compare_tsv("guilliams2022/13_annotation_of_record/r3_evidence.tsv")'''),
    *step(f"{G}/13b_aor_apply.py", 'run("13b_aor_apply.py", "rapids-sc", gpu=0)', "13b - annotation of record"),
    code('''compare_tsv("guilliams2022/13_annotation_of_record/composition_aor.tsv")
show("guilliams2022/13_annotation_of_record/figures/umap_aor.png")
show("guilliams2022/13_annotation_of_record/figures/dotplot_tableS1_aor.png")'''),
])

# ---------------------------------------------------------------- 07
save("07_samap_inputs.ipynb", [
    header("07 - SAMap inputs and the gene-map sanity check",
           """* SAMap inputs: mouse = annotation of record without 'unresolved', at most 4,000 cells per type (53,370); human = HLiCA v1.0,
  at most 1,200 per author_cell_type (50,400); var_names renamed to the SAMap gene names of the MMseqs2 maps (heavy notebook 11;
  record linked).
* Gene maps vs MGI curated 1:1 orthologs (Ensembl BioMart and HGNC HCOP could not be downloaded from the VM): result of record
  90.1 % (human->mouse) and 91.1 % (mouse->human) top hit, 97.6 % among the hits.""", "~10 min", "samap"),
    SETUP,
    code('''link("guilliams2022/13_annotation_of_record", "samap_maps/liver_hsmm/refs/mappings", "samap_maps/liver_hsmm/maps_mmseqs")
fresh("samap_run/input")'''),
    *step(f"{SR}/01_prep.py", 'run("01_prep.py", "samap")', "SAMap prep"),
    code('compare_tsv("samap_run/input/mm_subsample_counts.tsv")\ncompare_tsv("samap_run/input/hs_subsample_counts.tsv")'),
    code('fresh("samap_maps/liver_hsmm/refs/orthologs")\nlink("samap_maps/liver_hsmm/refs/orthologs/HOM_MouseHumanSequence.rpt")'),
    *step(f"{SM}/06_ortholog_check.py", 'run("06_ortholog_check.py", "samap")', "Ortholog check"),
    code('compare_tsv("samap_maps/liver_hsmm/refs/orthologs/ortholog_check.tsv", index_col=None)'),
])

# ---------------------------------------------------------------- 08
save("08_figures_tables.ipynb", [
    header("08 - Figures 1-3 and Table 1",
           """Fig 1 mouse atlas, Fig 2 HLiCA (authors' colour palette), Fig 3 SAMap (unintegrated row first), Table 1 SAMap matches
(kableExtra). Nature format: 183 x 170 mm, Arial 5-7 pt, audited with check_nature_fonts.py. The SAMap result itself comes from
heavy notebook 09 (record linked). The unintegrated UMAP is recomputed on GPU and may differ in layout from the record.""",
           "~15 min", "rapids-sc (figure data), samap (figures), masld_tables (R kableExtra), latex (tectonic), masld_atlas (audit)"),
    SETUP,
    code('''link("guilliams2022/13_annotation_of_record", "samap_run/input", "samap_run/out",
     "samap_maps/liver_hsmm/refs/orthologs/HOM_MouseHumanSequence.rpt")
fresh("figures", "figures/data", "figures/tables")'''),
    *step(f"{FG}/common.py", None, "Shared figure settings (palettes, lineage maps, mm layout helpers)"),
    *step(f"{FG}/fig_00_data.py", 'run("fig_00_data.py", "rapids-sc", gpu=0)', "Figure data: unintegrated baseline, dot plot statistics"),
    *step(f"{FG}/fig_00b_hlica_data.py", 'run("fig_00b_hlica_data.py", "samap")', "HLiCA figure data"),
    *step(f"{FG}/fig1_mouse_atlas.py", 'run("fig1_mouse_atlas.py", "samap")', "Figure 1 - mouse atlas"),
    *step(f"{FG}/fig2_hlica.py", 'run("fig2_hlica.py", "samap")', "Figure 2 - HLiCA"),
    *step(f"{FG}/fig3_samap.py", 'run("fig3_samap.py", "samap")', "Figure 3 - SAMap"),
    *step(f"{FG}/tab1_data.py", 'run("tab1_data.py", "samap")', "Table 1 data"),
    *step(f"{FG}/tab1_kable.R", 'run("tab1_kable.R", "masld_tables", cwd=rep("figures/tables"))', "Table 1 - kableExtra"),
    code('''import subprocess
print(subprocess.run([f"{ENVS}/latex/bin/tectonic", "-X", "compile", "table1_samap_matches.tex"], cwd=rep("figures/tables"),
                     capture_output=True, text=True).stderr[-400:])
compare_tsv("figures/tables/table1_samap_matches.tsv", index_col=None)'''),
    *step(f"{FG}/check_nature_fonts.py", 'run("check_nature_fonts.py", "masld_atlas", args=[rep(f"figures/{f}.pdf") for f in ("fig1_mouse_atlas", "fig2_hlica", "fig3_samap")])',
          "Nature format audit (all text 5-8 pt, page <= 183 x 170 mm)"),
    code('for f in ("fig1_mouse_atlas", "fig2_hlica", "fig3_samap"): show(f"figures/{f}.png", width=1000)'),
])

# ---------------------------------------------------------------- heavy
HEAVY_NOTE = """**Heavy notebook.** Nothing runs until you set `RUN_HEAVY = True` in the cell below. Long steps start in a detached
`screen` session with an exit sentinel (`WORK/logs/<name>.exit`); the wait cell polls the log. Re-run the light notebook named
below afterwards so that it uses these rebuilt outputs instead of the record's."""

save("09_heavy_samap_run.ipynb", [
    header("09 (heavy) - SAMap run, HLiCA x mouse", f"""sc-samap 3.0.1 (Project 23 WP2 recipe): SAMAP(filenames, f_maps), sm.run(), get_mapping_scores(n_top=0),
GenePairFinder(align_thr=0.1). Record: 62 min on 16 CPUs; outputs ~23 GB (pickled object 14 GB).

{HEAVY_NOTE} Afterwards: notebook 08.""", "~1 h, 16 CPU, ~25 GB disk", "samap"),
    SETUP, code("RUN_HEAVY = False   # set True to run"),
    code('''if RUN_HEAVY:   # nothing is touched unless the notebook is really run (review fix)
    guard(40)
    link("samap_run/input", "samap_maps/liver_hsmm/maps_mmseqs")
    fresh("samap_run/out")'''),
    *step(f"{SR}/02_samap.py", None),
    code('''launch("samap", [f"{PY['samap']} -I {WORK}/02_samap.py > {LOGS}/02_samap.py.log 2>&1"])
wait_for("samap", log=f"{LOGS}/02_samap.py.log")'''),
    code('''if RUN_HEAVY:
    import pandas as pd
    a = pd.read_csv(rec("samap_run/out/mapping_table.tsv"), sep="\\t", index_col=0); b = pd.read_csv(rep("samap_run/out/mapping_table.tsv"), sep="\\t", index_col=0)
    print("max |difference| of mapping scores vs record:", float((a - b.reindex_like(a)).abs().max().max()))
    # 02_samap.py logs a GenePairFinder failure instead of stopping; fail here if the gene pairs are missing (review fix)
    assert os.path.exists(rep("samap_run/out/gene_pairs_thr0.1.tsv")), "GenePairFinder failed: see the 02_samap log"'''),
])

save("10_heavy_doublets_decontx.ipynb", [
    header("10 (heavy) - scDblFinder and decontX per sample", f"""* 06a: split the raw object into one h5ad per sample on local disk (6.2 GB).
* 06b: scDblFinder (all cells, defaults, seed 12345) + default decontX (draft-QC cells, z = TotalVI Leiden r1), 46 samples, 12 in parallel.
* 06d: decontX strength ladder on 10 samples: default vs fixed delta (1,100), (1,300), (1,1000).
* 06e: final decontX delta (1,300), estimateDelta FALSE, for the 13 heavy in vivo cell samples.
R scDblFinder + celda through rpy2. Record: 06b ~20 min, 06d ~10 min, 06e ~25 min.

{HEAVY_NOTE} Afterwards: notebooks 04 onwards.""", "~1 h, 12-13 parallel R workers", "rapids-sc (06a), sc-preprocessing (06b, 06d, 06e)"),
    SETUP, code("RUN_HEAVY = False   # set True to run"),
    code('''SPLIT = os.path.expanduser("~/liver_cross_species/repro_tmp_split"); G6 = rep("guilliams2022/06_doublets_ambient")
if RUN_HEAVY:   # nothing is touched unless the notebook is really run (review fix)
    guard(20, os.path.expanduser("~")); os.makedirs(SPLIT, exist_ok=True)
    link("guilliams2022/01_raw", "guilliams2022/02_qc", "guilliams2022/05_annotation_draft", "guilliams2022/03_integration")
    fresh("guilliams2022/06_doublets_ambient/per_sample", "guilliams2022/06_doublets_ambient/ladder", "guilliams2022/06_doublets_ambient/final_dx")'''),
    *step(f"{G}/06a_split_samples.py", 'if RUN_HEAVY: run("06a_split_samples.py", "rapids-sc", args=[rep("guilliams2022"), SPLIT])'),
    *step(f"{G}/06b_doublets_decontx_sample.py", None),
    *step(f"{G}/06d_decontx_ladder.py", None),
    *step(f"{G}/06e_decontx_final_sample.py", None),
    code('''P = PY["sc-preprocessing"]; T = "export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2"
LADDER = "CS140 CS139 CS93 CS144 CS89 CS115 ABU19 ABU20 ABU21 CS98"
TIER1 = "CS114 CS137 CS138 CS139 CS144 CS48 CS52 CS88 CS89 CS92 CS93 CS96 CS97"
launch("dbl_dx", [T, f"mkdir -p {LOGS}/06b {LOGS}/06d {LOGS}/06e",
    f"ls -S {SPLIT}/*.h5ad | xargs -n1 basename | sed 's/\\\\.h5ad$//' | xargs -P 12 -I{{}} sh -c '{P} -I {WORK}/06b_doublets_decontx_sample.py {SPLIT} {G6}/per_sample {{}} > {LOGS}/06b/{{}}.log 2>&1'",
    f"echo {LADDER} | tr ' ' '\\\\n' | xargs -P 10 -I{{}} sh -c '{P} -I {WORK}/06d_decontx_ladder.py {SPLIT} {G6}/ladder {{}} > {LOGS}/06d/{{}}.log 2>&1'",
    f"echo {TIER1} | tr ' ' '\\\\n' | xargs -P 13 -I{{}} sh -c '{P} -I {WORK}/06e_decontx_final_sample.py {SPLIT} {G6}/final_dx {{}} > {LOGS}/06e/{{}}.log 2>&1'",
    f"[ $(ls {G6}/per_sample/*.cells.parquet | wc -l) -eq 46 ] && [ $(ls {G6}/final_dx/*.dx300.npz | wc -l) -eq 13 ]"])
wait_for("dbl_dx")'''),
    md("Delete the split folder (`~/liver_cross_species/repro_tmp_split`, 6.2 GB) when notebook 04 has been re-run."),
])

save("11_heavy_gene_maps_mmseqs2.ipynb", [
    header("11 (heavy) - Human-mouse gene maps for SAMap (MMseqs2)", f"""Copies of the Project 23 WP2 scripts (https://github.com/fbastian/BH2026-23-single-cell, commit in UPSTREAM_COMMIT.txt;
our changes in 02_vs_upstream.diff / 03_vs_upstream.diff): GENCODE v48 (HLiCA) and vM18 (= Ensembl 93, Guilliams mm10-3.0.0)
transcripts, one transcript per gene, SAMap names (Symbol or Symbol__EnsemblID), MMseqs2 18.8cc5c translated search both
directions (-s 7.5 -e 1e-6 --min-length 30 ...). Record: ~45 min per direction on 20 threads.

**Disk: the search needs ~170 GB of temporary space under ROOT/tmp (on 11 Nov the vault dropped to 15 GB free).** The guard
below refuses to start with less than 200 GB free.

{HEAVY_NOTE} Afterwards: notebooks 07 and 09.""", "~1.5 h + downloads, 20 threads, >= 200 GB temporary disk",
           "rapids-sc (+ mmseqs2 binary from env mmseqs2)"),
    SETUP, code("RUN_HEAVY = False   # set True to run"),
    code('''if RUN_HEAVY:   # nothing is touched unless the notebook is really run (review fix)
    guard(200)
    link("guilliams2022/01_raw")
    fresh("samap_maps/liver_hsmm/refs/raw", "samap_maps/liver_hsmm/refs/samap", "samap_maps/liver_hsmm/refs/mappings", "samap_maps/liver_hsmm/maps_mmseqs")'''),
    *step(f"{SM}/01_download.sh", 'if RUN_HEAVY: run("01_download.sh", "rapids-sc")'),
    code('''if RUN_HEAVY:   # 01_download.sh prints BAD but exits 0 on a broken file; check against the record checksums (review fix)
    rec_sums = dict(l.split()[::-1] for l in open(rec("samap_maps/liver_hsmm/refs/raw/SHA256SUMS")))
    new_sums = dict(l.split()[::-1] for l in open(rep("samap_maps/liver_hsmm/refs/raw/SHA256SUMS")))
    assert new_sums == rec_sums, f"downloads differ from the record: {new_sums} vs {rec_sums}"; print("downloads match the record checksums")'''),
    *step(f"{SM}/02_build_samap_fastas_liver.py", 'if RUN_HEAVY: run("02_build_samap_fastas_liver.py", "rapids-sc")'),
    *step(f"{SM}/03_build_mmseqs_maps_liver.py", None),
    code('''launch("mmseqs", [f"PATH={ENVS}/mmseqs2/bin:$PATH {PY['rapids-sc']} -I {WORK}/03_build_mmseqs_maps_liver.py > {LOGS}/03_mmseqs.log 2>&1"])
wait_for("mmseqs", log=f"{LOGS}/03_mmseqs.log", every=300)'''),
])

save("12_heavy_totalvi.ipynb", [
    header("12 (heavy) - TotalVI, three times", f"""The authors' integration (4c_TotalVI_CiteSeq+non-CiteSeq.ipynb): 4000 HVG (seurat, batch-aware by sample), TOTALVI n_latent 20,
encoder covariates, lr 1e-3, 400 epochs, train 0.9, batch = sample (scvi-tools 1.5; batch_size 1024, 8 loader workers).
* 03 + 03p: draft (02_qc, 328,041 cells). Record 114 min.
* 08 + 08p: clean counts (07_clean, 321,059). Record 78 min.
* 11 + 11p: strict cells (10_strict, 234,552). Record 112 min.
Post steps (neighbours, UMAP, Leiden) run on GPU with rapids-singlecell. Leiden numbers can differ from the record; the hand
decisions in notebooks 04-06 name record clusters, so compare before re-using them (ARI cell below).

{HEAVY_NOTE} Afterwards: notebooks 03 onwards. Each training run needs its input from the previous light notebooks (or the
record): 03 <- 02, 08 <- 04/05 (07_clean), 11 <- 05 (10_strict).""", "~5 h on one T4 GPU", "masld_integration (training), rapids-sc (post)"),
    SETUP, code("RUN_HEAVY = False   # set True to run"),
    code('''if RUN_HEAVY:   # nothing is touched unless the notebook is really run (review fix)
    guard(30)
    link("guilliams2022/02_qc", "guilliams2022/07_clean", "guilliams2022/10_strict")
    fresh("guilliams2022/03_integration", "guilliams2022/08_totalvi_clean", "guilliams2022/11_totalvi_strict")'''),
    *step(f"{G}/03_integrate.py", None, "Draft TotalVI"), *step(f"{G}/03p_post_gpu.py", None),
    *step(f"{G}/08_totalvi_clean.py", None, "TotalVI on the clean counts"), *step(f"{G}/08p_post_gpu.py", None),
    *step(f"{G}/11_totalvi_strict.py", None, "TotalVI on the strict cells"), *step(f"{G}/11p_post_gpu.py", None),
    code('''MI, RS, O = PY["masld_integration"], PY["rapids-sc"], rep("guilliams2022")
launch("totalvi", ["export CUDA_VISIBLE_DEVICES=0", "set -e",
    f"{MI} -I {WORK}/03_integrate.py {O} totalvi > {LOGS}/03_integrate.py.log 2>&1",
    f"{RS} -I {WORK}/03p_post_gpu.py {O} > {LOGS}/03p_post_gpu.py.log 2>&1",
    f"{MI} -I {WORK}/08_totalvi_clean.py > {LOGS}/08_totalvi_clean.py.log 2>&1",
    f"{RS} -I {WORK}/08p_post_gpu.py > {LOGS}/08p_post_gpu.py.log 2>&1",
    f"{MI} -I {WORK}/11_totalvi_strict.py > {LOGS}/11_totalvi_strict.py.log 2>&1",
    f"{RS} -I {WORK}/11p_post_gpu.py > {LOGS}/11p_post_gpu.py.log 2>&1"])
wait_for("totalvi", every=600)'''),
    code('''if RUN_HEAVY:
    import anndata as ad
    from sklearn.metrics import adjusted_rand_score as ari
    for f, k in [("03_integration/embeddings.h5ad", "leiden_totalvi_r1"), ("08_totalvi_clean/embeddings.h5ad", "leiden_clean_r3"),
                 ("11_totalvi_strict/embeddings.h5ad", "leiden_clean_r3")]:
        a = ad.read_h5ad(rec(f"guilliams2022/{f}"), backed="r").obs[k]; b = ad.read_h5ad(rep(f"guilliams2022/{f}"), backed="r").obs[k]
        print(f, k, "ARI vs record", round(ari(a.astype(str), b.reindex(a.index).astype(str)), 3))'''),
])
