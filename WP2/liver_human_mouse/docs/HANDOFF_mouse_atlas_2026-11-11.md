# Handoff: Guilliams 2022 Fig 1B rebuild + human–mouse liver maps (BioHackathon 2026, Project 23)

Written 2026-11-11 13:00 CET. Times in CET (VM clock is UTC, +1 h).

## Goal

Cross-species (human + mouse) liver integration as groundwork for MASLD cross-species work.
Step 1: rebuild the healthy mouse liver atlas of Guilliams et al., Cell 2022 (Fig 1B, 185,894 cells/nuclei,
TotalVI) from GEO, annotate it ourselves, then integrate with human liver (HLiCA) using SAMap gene maps built
with the Project 23 WP2 scripts (https://github.com/fbastian/BH2026-23-single-cell).

## Where things are

| What | Path |
|---|---|
| Code (this project) | `~/liver_cross_species/guilliams2022/code/` (numbered steps + `run_*.sh` screen launchers) |
| Logs + exit sentinels | `~/liver_cross_species/guilliams2022/logs/` (`<step>.exit` = exit code) |
| Outputs (vault: `/` is below its 200 GB guard) | `/mnt/masld_vault/cross_species/guilliams2022/` |
| GEO single-cell data (all species) | `/mnt/masld_vault/geo/GSE192740/` (RAW.tar, sample table, human Seurat rds) |
| GEO spatial (downloaded, not used) | `/mnt/masld_vault/geo/GSE192741/` |
| Paper supplements (Tables S1–S10) | `/mnt/masld_vault/geo/GSE192740/paper_supplement_cell2022/` |
| Paper PDF | `~/MASLD_atlas/docs/atlas_papers/guilliams_2022_liver_macrophage_niches.pdf` |
| Author code | `/mnt/masld_vault/geo/github-guilliottslab-scripts_GuilliamsEtAll_Cell2022/` |
| HLiCA full object (human) | `/mnt/masld_vault/geo/cxg-059202e1-1f1b-483f-9151-f3a25a380c39/HLiCA_human_liver_all_cells.h5ad` |
| Hackathon repo clone | `~/liver_cross_species/BH2026-23-single-cell/` |

## Data facts (verified)

- Fig 1B = 46 GEO samples: mouse rows titled `Healthy*`, `PBS Oral Gavage*`, `Standard diet fed (12 weeks)*`
  → 15 ex vivo / 19 in vivo / 12 nuclei, 18 CITE-seq. Matches the Fig 1A legend. 389,056 cells before QC.
- Cell counts per sample = GEO table for all 46. ADT rows per CITE sample = GEO "number of added ABs"
  (107/114/115/169); 173 antibodies after the authors' name translation (`4f_TranslateABs.txt`).
- CITE-seq samples are dense txt.gz (genes then ADT rows, R `make.unique` suffixes); others 10x h5, genes only.
- Mouse reference: 10x `mm10-3.0.0.premrna` (GRCm38, Ensembl 93), 31,053 genes; `var_names` = Ensembl IDs,
  symbols in `var['symbol']` (71 genes share 36 duplicated symbols).
- No mouse object or per-cell labels on GEO or in the supplements. Labels (`annot_mouseStStAll.csv`) were on
  livercellatlas.org, offline since a PHP security issue; Wayback has only `annot_digestExVivo.csv`.
- Table S1 (mmc1) = DE genes + proteins for exactly the 17 Fig 1B cell types → our annotation reference.
- GEO human rds = the paper's Fig 4B (167,598 cells, `annot` 18 types, RNA counts + denoised ADT).
  Double-gzipped; Seurat 3.1.2, `UpdateSeuratObject` fails validation → read slots directly.
- HLiCA: CELLxGENE schema 7.1.0, 524,699 cells × 32,624 genes, raw counts in `raw.X`, Ensembl IDs as
  `var_names`, symbols in `var['feature_name']` (20 duplicated).

## Pipeline status (13:00)

| Step | Script | Status | Output |
|---|---|---|---|
| 01 load 46 samples | `01_load_raw.py` | done | `01_raw/guilliams2022_fig1b_raw.h5ad` (389,056 × 31,053, `obsm['adt']`, `obsm['adt_measured']`) |
| 02 draft QC (authors' rule) | `02_qc.py` | done (GPU Scrublet) | `02_qc/` – keeps 328,041 (paper 185,894) |
| 02b pre-QC L0 + metrics | `02b_preqc_L0_metrics.py` | done | `02_qc_v1/cells_preqc.parquet` |
| 02c QC v1 thresholds + ridges | `02c_thresholds_ridges.py` | done, **not applied** | `02_qc_v1/figures/` (12 ridges), would remove 10.8 % |
| 03h unintegrated + Harmony | `03h_harmony_gpu.py` | done (4.5 min) | `03_integration/embeddings_harmony.h5ad` |
| 03 TotalVI (draft) | `03_integrate.py` → `03p_post_gpu.py` → `04_umap_draft_panel.py` | running, screen `g22_totalvi`, epoch 30/400, 15.9 s/epoch | ETA ~14:35 train, ~14:45 panel |
| Draft panel | `04_umap_draft_panel.py` | Harmony rows done; TotalVI row added automatically | `03_integration/figures/umap_draft_panel.pdf` |

## Known issues / deviations

1. **No ambient correction.** Authors used FastCAR (needs unfiltered droplets; GEO has filtered only). Alb is
   high in nearly every cluster on the draft UMAP. Candidate fix: decontX (as MASLD v1).
2. **QC v1 lineage labels need fixing before locking.** Marker-score L0 is inflated by ambient hepatocyte RNA:
   epithelial strata in CS53/CS62/CS115/CS136/CS138 are bimodal and the 5-MAD upper bound would cut real
   hepatocytes; only 730 of 37,138 nuclei are called immune. Fix: ambient-aware L0, drop epithelial upper bound.
3. **Doublets:** Scrublet makes almost no calls (as in MASLD, where it made zero in all Guilliams human
   libraries). scDblFinder per sample is still to run; MASLD rule = remove only on Scrublet + scDblFinder consensus.
4. **TotalVI training settings:** batch 1024 + 8 loader workers (authors: scvi default 256). Draft only.
   First run at batch 256 was stopped (log `03_totalvi_bs256_stopped.log`).
5. Draft QC keeps far more cells than the paper (328k vs 186k): authors also removed contaminating /
   proliferating clusters by hand per sample.

## Gene maps for SAMap (WP2)

- Their refs: human Ensembl 93 (GRCh38), mouse GENCODE vM12 (= Ens87, GRCm38), zebrafish Ens116.
  Assemblies match ours; annotations do not (mouse Ens93 here; HLiCA newer GENCODE).
- Their maps (`maps/hsmm`, `refs/mappings/*.gene_names.tsv`) are on the project OneDrive, not GitHub, and cover
  only the genes of their muscle datasets. Requested on Slack 2026-11-11.
- Plan: run their `02_build_samap_fastas.py` + `03_build_mmseqs_maps.py` (same parameters) on our liver gene
  sets with matching cDNA; rename `var_names` to `Symbol` / `Symbol__EnsemblID`; compare with their `hsmm`.

## Gotchas learned today

- Never stop a screen job with `pgrep -f <script> | kill`: it also matches the SCREEN process and kills the child.
- rsc Scrublet on all 389k cells at once OOMs a T4; per sample it takes < 3 min in total.
- rsc scale + PCA of 328k × 4000 OOMs a T4 → scale/PCA on CPU, Harmony/neighbours/UMAP on GPU.
- Scrublet in `masld_integration` needs scikit-image (missing) → QC steps run in `rapids-sc`.
- Disk: `/` at ~190 GB free (below the 200 GB guard); vault 166 GB free. Slim-archiving old scAtlasTb test runs
  would free ~300 GB (waiting for the user's go).

## Update 13:10 – gene maps started (steps 9–10)

- Code: `~/liver_cross_species/samap_maps/code/` = copies of WP2 scripts 02 + 03 (upstream commit in
  `UPSTREAM_COMMIT.txt`, minimal diffs in `02_vs_upstream.diff`, `03_vs_upstream.diff`).
- References: human GENCODE v48 (HLiCA, CxG schema 7.1.0), mouse GENCODE vM18 (= Ensembl 93 = 10x mm10-3.0.0),
  in `/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/refs/raw` (+ SHA256SUMS).
- Step 9 done: SAMap names + one-transcript-per-gene FASTAs; coverage human 32,624/32,624, mouse 30,976/31,053
  (77 missing: clone-based AC… genes and PAR genes). Tables in `refs/mappings/{hs,mm}.*.tsv`.
- Step 10 running: MMseqs2 v18, their parameters, screen `hsmm_mmseqs`, output `maps_mmseqs/hsmm/{hs_to_mm,mm_to_hs}.txt`.
