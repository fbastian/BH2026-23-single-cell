# Handoff: human–mouse liver gene maps for SAMap (BioHackathon 2026, Project 23)

Written 2026-11-11 14:45 CET. Times in CET (VM clock is UTC, +1 h).
Companion handoff for the mouse atlas: `~/liver_cross_species/guilliams2022/docs/HANDOFF_2026-11-11.md`.

## Goal

SAMap needs a gene homology graph between the two species: a list of human gene ↔ mouse gene pairs with alignment
scores. We built that graph for the two liver objects we want to integrate:

- **Human:** HLiCA v1.0 (CELLxGENE, 524,699 cells × 32,624 genes)
- **Mouse:** Guilliams et al., Cell 2022 healthy atlas (Fig 1B samples, 31,053 genes)

We used the Project 23 WP2 scripts and settings (https://github.com/fbastian/BH2026-23-single-cell) on **our** gene
sets. Their own maps only cover the genes in their muscle datasets.

## Status

| Step | What | Status |
|---|---|---|
| 1 | Download matching transcript FASTAs | done |
| 2 | One transcript per gene + SAMap gene names (`02_build_samap_fastas_liver.py`) | done |
| 3 | MMseqs2 translated search both directions (`03_build_mmseqs_maps_liver.py`) | **done ~14:30 CET**, exit 0 |
| 4 | Rename `var_names` of both objects to the SAMap names | not started |
| 5 | Compare with the team's `hsmm` map (OneDrive) | waiting: team has not shared `maps/hsmm` + `refs/mappings` yet |
| 6 | Sanity check against a curated ortholog table (Ensembl BioMart / HGNC HCOP) | not started |
| 7 | SAMap run HLiCA × Guilliams | not started (Thu) |

## Where things are

| What | Path |
|---|---|
| Code | `~/liver_cross_species/samap_maps/code/` |
| Logs + exit sentinel | `~/liver_cross_species/samap_maps/logs/` (`03_mmseqs.log`, `03.exit`) |
| Data root (vault) | `/mnt/masld_vault/cross_species/samap_maps/liver_hsmm/` |
| Raw references + checksums | `refs/raw/` (`gencode.v48.transcripts.fa.gz`, `gencode.vM18.transcripts.fa.gz`, `SHA256SUMS`) |
| SAMap FASTAs | `refs/samap/human.GRCh38.GENCODEv48.samap.fa`, `refs/samap/mouse.GRCm38.GENCODEvM18.samap.fa` |
| Gene name tables | `refs/mappings/{hs,mm}.gene_names.tsv` (`ensembl_gene_id`, `gene_symbol`, `samap_gene`, `symbol_ambiguous`) |
| Chosen transcripts / missing genes | `refs/mappings/{hs,mm}.selected_transcripts.tsv`, `{hs,mm}.missing_from_transcriptome.tsv` |
| **Gene maps (SAMap input)** | `maps_mmseqs/hsmm/hs_to_mm.txt`, `maps_mmseqs/hsmm/mm_to_hs.txt`, `mmseqs_version.txt` |
| Hackathon repo clone | `~/liver_cross_species/BH2026-23-single-cell/` |

## Scripts: copied from the team's repo, not rewritten

- Scripts 02 and 03 are copies of `WP2/scripts/02_build_samap_fastas.py` and `03_build_mmseqs_maps.py`. The
  upstream commit is in `code/UPSTREAM_COMMIT.txt` (`85d99fba…`), and every change is in `02_vs_upstream.diff` and
  `03_vs_upstream.diff`.
- **Changed:**
  - Paths and the species pair (human + mouse only, no zebrafish).
  - Our h5ads and their gene fields: HLiCA `var['feature_name']`, Guilliams `var['symbol']`.
  - The GENCODE FASTA parser also accepts human `ENSG` IDs.
  - `THREADS = 20`, so TotalVI keeps CPU on the shared VM.
- **Unchanged:** gene naming, transcript choice and all search parameters.
- `01_download.sh` and `run_03.sh` are ours.

## Reference choices

| | Team (muscle) | Ours (liver) | Why |
|---|---|---|---|
| Human | Ensembl 93 cDNA | GENCODE v48 transcripts | HLiCA = CELLxGENE schema 7.1.0 → GENCODE v48 gene IDs |
| Mouse | GENCODE vM12 (= Ens87) | GENCODE vM18 (= Ensembl 93) | Guilliams used 10x `mm10-3.0.0.premrna` = Ensembl 93 |

The assemblies match the team's (GRCh38, GRCm38); the annotation releases differ, so their maps cannot be reused
as they are.

## Step 2 result: SAMap gene names

- **Naming rule (team's):** `samap_gene` = symbol. When a symbol is shared by several gene IDs, it is
  `Symbol__EnsemblID`.
- **Human:** 32,624 / 32,624 genes have a transcript (100 %); 20 ambiguous symbols.
- **Mouse:** 30,976 / 31,053 (99.75 %); 71 ambiguous. The 77 missing are clone-based `AC…` genes and PAR genes
  that are absent from vM18 transcripts. Listed in `mm.missing_from_transcriptome.tsv`.

## Step 3 result: MMseqs2 maps

- **Tool:** MMseqs2 `18.8cc5c`, conda env `mmseqs2`.
- **Search:** `mmseqs search --search-type 2` (translated nucleotide vs nucleotide), with:
  - `-s 7.5 -e 1e-6 --min-length 30 --orf-start-mode 1`
  - `--max-seqs 10000 --max-accept 1000000 --max-rejected 1000000 -a`
- **Post-processing:** best alignment kept per gene pair. Output columns (no header): query, target, pident, alnlen,
  mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits.
- **Runtime:** about 45 min per direction (prefilter ~15 min + align ~27 min), ~1.5 h in total on 20 threads.

| Map | Pairs | Query genes with ≥ 1 hit | Median hits per gene | Top hit has the same symbol |
|---|---|---|---|---|
| `hs_to_mm.txt` | 1,029,058 | 23,029 / 32,624 | 4 | 58.4 % |
| `mm_to_hs.txt` | 1,011,410 | 22,560 / 31,053 | 6 | 58.7 % |

**Liver marker spot check** (human → best mouse hit, identity %, bits):

| Human | Mouse | Identity | Bits |
|---|---|---|---|
| ALB | Alb | 72.4 % | 964 |
| CLEC4F | Clec4f | 49.1 % | 484 |
| PTPRC | Ptprc | 73.0 % | 772 |
| STAB2 | Stab2 | 78.1 % | 4426 |
| CYP2E1 | Cyp2e1 | 78.7 % | 475 |

The "same symbol" rate is not an accuracy figure. The gene sets include lncRNAs and other non-coding genes, whose
best hit is often a paralog or a `Gm`/`Rik` gene. Check it properly against a curated ortholog table restricted to
protein-coding genes (step 6) before we quote it.

## Next steps

1. **Rename objects (step 4):** set `var_names` of HLiCA and the Guilliams object to `samap_gene` using
   `{hs,mm}.gene_names.tsv`, and drop the 77 mouse genes without a transcript. SAMap matches map names to
   `var_names` exactly.
2. **Ortholog sanity check (step 6):** compare the top hits with Ensembl BioMart 1:1 orthologs (protein-coding).
3. **Compare with the team (step 5):** when they share `maps/hsmm`, compare pair overlap and bit scores on the
   genes both maps cover.
4. **SAMap (step 7):** HLiCA × Guilliams, first on the draft-QC mouse object; rerun on the locked QC v1 object.
   Subsample HLiCA, because SAMap on 525k + 328k cells is heavy.

## Gotchas

- **Disk:** the translated search needs about **170 GB of temporary space** (~82 GB per direction), written under
  `ROOT/tmp`, which is on the vault. On 11 Nov the vault dropped to 15 GB free before the script cleaned up after
  itself. Before any rerun (e.g. adding zebrafish), check free space, or point `ROOT/tmp` somewhere with ≥ 200 GB.
- The script refuses to run if `maps_mmseqs/hsmm` already exists. Use a new `OUT` for a new run.
- MMseqs2 accession parsing truncates custom names. The team's script works around it by exporting full headers;
  keep that if editing.
- Gene-map IDs are SAMap names (`Symbol` / `Symbol__ENS…`), not Ensembl IDs. Our h5ads currently have Ensembl IDs as
  `var_names`.
