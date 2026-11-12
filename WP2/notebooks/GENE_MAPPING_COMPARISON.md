# Gene-mapping method comparison

**BioHackathon Europe 2026 — Project 23, WP2**

**Contributor:** Nada Orsulics

## Overview

This notebook compares precomputed SAMap MappingT matrices
for three gene-mapping strategies:

- BLAST (map1)
- MMseqs2 (map2)
- BLAST + OrthoFinder (map3)

The preliminary analysis uses mouse and zebrafish cell types.

## Analysis

The notebook includes:

- Pairwise Spearman correlations and mean absolute errors
- MappingT heatmaps
- Top-1 cross-species matches
- Top-3 match stability
- Pairwise MappingT difference heatmaps

These comparisons measure agreement between methods,
not biological mapping accuracy.

## Input files

The following files are required:

- MappingT_mm_dr_annot.csv
- MappingT_mm_dr_annot_map2.csv
- MappingT_mm_dr_annot_orthofinder.csv

Input matrices are not included in this contribution.
They must be obtained from the WP2 project data.

## Running the notebook

Required Python packages:

- numpy
- pandas
- scipy
- matplotlib
- nbformat
- jupyter

From the repository root, run:

```bash
SAMAP_RESULTS_DIR="/path/to/MappingT/files" \
jupyter nbconvert --to notebook --execute \
WP2/notebooks/compare_gene_mapping_methods.ipynb \
--output compare_gene_mapping_methods_executed.ipynb
```

## Limitations

Results are preliminary and based on precomputed MappingT
matrices. The two cross-species mapping directions are
related observations, not independent biological replicates.

The notebook can be rerun with updated matrices after
full-dataset analyses, including human data if the
input handling is extended accordingly.
