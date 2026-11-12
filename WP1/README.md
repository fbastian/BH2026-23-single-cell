# WP1: Biodata resource integration

Pipeline to build the dataset metadata used for the species x factor table, online at **https://fbastian.github.io/BH2026-23-single-cell/**. The commands of each step are run from this folder (`WP1/`), in the conda environment described in [SETUP.md](SETUP.md). Input and intermediate files live in `data/` (git-ignored); generated files to share, such as the term remapping TSV, go in `output/`, and the files loaded by the explore page in `explore/ui/public/`.

To run steps 1 to 6 in one go, in the conda environment (the script can be started from any folder, here the repository root):

```
bash WP1/scripts/run_pipeline.sh                  # downloads the ontologies again
bash WP1/scripts/run_pipeline.sh --skip-download  # uses the ontologies already in data/
```

[scripts/run_pipeline.sh](scripts/run_pipeline.sh) runs the commands below and stops at the first error. The dataset exports (`data/datasets_CxG.json`, `data/datasets_scFAIR.json`, `data/datasets_fenix.json`) must already be in `data/`. Set `PYTHON` to use another interpreter (e.g. `PYTHON=/path/to/envs/bh23/bin/python`). A new ontology release can change the remappings: check the warnings and the `unresolved` rows of step 3.

## 1. Download the ontologies

From the `data` folder:

```
curl -sLO http://purl.obolibrary.org/obo/uberon/composite-metazoan.obo
curl -sLO http://purl.obolibrary.org/obo/uberon/collected-metazoan.obo
curl -sLO http://purl.obolibrary.org/obo/cl.obo
curl -sLO http://purl.obolibrary.org/obo/zfa.obo
```

`composite-metazoan.obo` is the reference ontology. The other three are only used to find the closest ancestors of terms that composite does not know (step 3).

## 2. Merge the datasets of the different sources

[scripts/merge_datasets.py](scripts/merge_datasets.py) merges the JSON exports of several sources into one deduplicated JSON.

```
python scripts/merge_datasets.py \
  --input cxg=data/datasets_CxG.json \
  --input scfair=data/datasets_scFAIR.json \
  --input fenix=data/datasets_fenix.json \
  --out output/datasets_merged.json \
  --report output/merged_datasets.tsv
```

- `--input FORMAT=PATH` is repeated once per source. `FORMAT` selects how the file is read: `cxg` (CELLxGENE curation API export), `scfair` (scFAIR portal export) or `fenix` (Single Cell Expression Atlas export). To add a source, write a reader function in the script and register it in `READERS`.
- **The order of the `--input` options is the priority**: the first source is preferred for scalar fields. Put the freshest source first.
- Deduplication:
  - the same dataset of the same source (same source and `source_reference_id`) is always merged;
  - datasets of different resources with the same DOI, the same organisms and the same cell count are merged too (only if unambiguous), and flagged with `"merged": true`. `--report` lists these merges.
- Merged records keep the most complete information: lists (tissue, cell type, ...) are the union of the sources. A record of an older version of a dataset (different `dataset_version_id`) only fills empty fields. scFAIR copies of CELLxGENE datasets have no version id: it is read from the name of their H5AD file (`<dataset_version_id>.h5ad`).
- Each output record has a `dataset_id` (`<source>:<source id>`) and a `sources` list with the provenance of every merged record.
- DOIs are stored without URL encoding (Single Cell Portal gives e.g. `10.1016%2Fj.cell.2020.04.035`) and without `doi.org` prefix, and compared case-insensitively.

## 3. Generate the term remapping file

[scripts/make_term_remappings.py](scripts/make_term_remappings.py) lists the annotations (tissue, cell type, developmental stage) of the merged datasets that must be remapped to terms of the composite-metazoan Uberon, and writes them to a TSV.

```
python scripts/make_term_remappings.py \
  --composite data/composite-metazoan.obo \
  --datasets output/datasets_merged.json \
  --out output/term_remappings.tsv \
  --extra-obo data/collected-metazoan.obo data/cl.obo data/zfa.obo \
  --manual term_remappings_manual.tsv
```

`--composite`, `--datasets` and `--out` are mandatory; the script downloads nothing. The remapping rules, in this order:

1. `unknown`, `na` or an empty annotation is remapped to the root of the ontology (tissue `UBERON:0001062`, cell type `CL:0000000`, developmental stage `UBERON:0000105`).
2. Cell lines (`CVCL_*`) are excluded.
3. A term that exists in composite is kept as is (no row). Stubs do not count: composite keeps an empty stanza (no name, no parents) for terms merged into another term, e.g. zebrafish `ZFA:0009248` merged into `CL:0000540` neuron; rule 5 remaps them.
4. An obsolete term with a replacement (`replaced_by`) is remapped to it (e.g. `CL:0000089` to `CL:0000016`).
5. A term that composite only knows as an xref of another term is remapped to that term (e.g. `HsapDv:0010000` to `UBERON:0000066`).
6. Any other term is remapped to its closest (most specific) ancestor present in composite and in the domain of the factor (see below), using the ontologies given with `--extra-obo`. Without `--extra-obo`, these terms get the action `unresolved`.

Manual remappings go in [term_remappings_manual.tsv](term_remappings_manual.tsv) (columns `factor`, `source_term_id`, `action` [`remap` or `exclude`], `target_term_id`, `comment`), given with `--manual`. They replace whatever the rules above produce for that term; a term mapped to several targets has one row per target. It currently holds the obsolete CL terms without replacement, the zebrafish segmentation and hatching stages, and `UBERON:0000104` (life cycle, used as a stage).

The output TSV has one row per remapping (columns: `factor`, `source_term_id`, `source_label`, `action`, `target_term_id`, `target_label`, `via_term_id`, `via_label`, `method`, `n_datasets`, `organisms`, `comment`). Rows with the action `unresolved` must be fixed in the manual file, not in the output TSV, which is overwritten each time the script runs.

When reading the ontologies, ancestors follow `is_a` and `part_of` (the FlyBase relation "substage of", `FBdv:00018001`, counts as `part_of`). Taxon-conditional relations (`{gci_relation="in_taxon", ...}`, e.g. embryo stage `is_a` "Brugia life stage") are ignored.

Because `part_of` crosses between kinds of entities, ancestors are restricted to the domain of each factor (this applies to rule 6 and to the ontology roll-up):

| Factor | Domain |
|---|---|
| cell type | `CL:0000000` (cell) and its descendants |
| tissue | `UBERON:0001062` (anatomical entity) and its descendants, except cell types (descendants of `CL:0000000`) |
| developmental stage | `UBERON:0000105` (life cycle stage) and its descendants |

Terms of BFO, NCBITaxon, COB and GO are never used. Without these restrictions, cell types would roll up into organs (e.g. "decidual natural killer cell, human" into decidua), and organs into taxa (Metazoa, through "multicellular organism"). A cell type annotated as tissue (e.g. a primary cell culture) stays a tissue column. It rolls up only through `part_of` into anatomy (e.g. tracheal epithelial cell to trachea), and it is remapped within the cell type domain if needed. The script warns about targets outside their domain, including xref and manual targets.

## 4. Expand the terms to their ancestors

[scripts/make_term_ancestors.py](scripts/make_term_ancestors.py) applies the remappings to the dataset annotations, then lists, for each factor, every annotated term and its ancestors in composite, following `is_a` and `part_of` within the domain of the factor (table above). The datasets JSON is not modified: the expansion is a separate file.

```
python scripts/make_term_ancestors.py \
  --composite data/composite-metazoan.obo \
  --datasets output/datasets_merged.json \
  --remappings output/term_remappings.tsv \
  --out output/term_ancestors.tsv
```

All options are mandatory. Use the same datasets file as in step 3. The TSV has one row per term and factor, for annotated terms and for their ancestors:

| Column | Content |
|---|---|
| `factor` | `tissue`, `cell_type` or `development_stage` |
| `term_id`, `term_label` | term of composite and its name |
| `annotated` | `yes` if a dataset is annotated with the term (after remapping), `no` if it is only an ancestor |
| `parent_ids` | closest ancestors within the domain (`;`-separated): the hierarchy, e.g. for drill-down |
| `ancestor_ids` | all ancestors within the domain (`;`-separated), the term excluded: a dataset annotated with the term counts for all of them in the roll-up |

The script warns about annotated terms that are not in composite (e.g. a remapping TSV generated from other datasets) or outside their factor's domain.

## 5. Tag the uninformative terms

[scripts/make_term_subsets.py](scripts/make_term_subsets.py) lists the terms of the ancestor TSV (step 4) that are in given subsets of composite (`subset:` lines of the OBO file). The explore page hides them in the filter trees and in the roll-up columns of the pivot, and shows their more specific terms instead.

```
python scripts/make_term_subsets.py \
  --composite data/composite-metazoan.obo \
  --ancestors output/term_ancestors.tsv \
  --subsets non_informative upper_level ubprop:upper_level common_anatomy \
  --manual term_subsets_manual.tsv \
  --out output/term_subsets.tsv
```

`--composite`, `--ancestors`, `--subsets` and `--out` are mandatory. The subsets:

- `non_informative` ("abstract class brought in to group ontology classes but not informative", e.g. regional part of brain);
- `upper_level` ("abstract upper-level terms not directly useful for analysis", e.g. material anatomical entity, organ);
- `ubprop:upper_level`, the same tag on CL terms (cell, germ line cell, stuff accumulating cell);
- `common_anatomy` ("terms applicable across life, not restricted to Metazoa"), mostly the same terms as `upper_level`, plus generic ones such as anatomical conduit, disconnected anatomical group, anatomical collection.

Hand-curated terms to hide go in [term_subsets_manual.tsv](term_subsets_manual.tsv) (columns `factor`, `term_id`, `comment`), given with `--manual`; they get the subset `manual`. A term that is not in the ancestor TSV is reported and ignored. The file is empty for now.

The output TSV has the columns `factor`, `term_id`, `term_label`, `annotated` (as in step 4) and `subsets`. Terms annotated to datasets can be in these subsets too (e.g. cortex, hippocampal field, multicellular anatomical structure, and the roots, which stand for missing annotations): they stay visible as raw terms in the pivot.

## 6. Export the data for the explore page

[scripts/export_explore.py](scripts/export_explore.py) writes the two Parquet files loaded by the explore page ([explore/](explore/)): the datasets with their remapped annotations, and the terms with their hierarchy.

```
python scripts/export_explore.py \
  --datasets output/datasets_merged.json \
  --remappings output/term_remappings.tsv \
  --ancestors output/term_ancestors.tsv \
  --subsets output/term_subsets.tsv \
  --out-datasets explore/ui/public/datasets.parquet \
  --out-terms explore/ui/public/terms.parquet
```

All options are mandatory. Use the files of steps 2 to 5. What the export does:

- Tissue, cell type and developmental stage are remapped with the TSV of step 3, like in step 4: an empty annotation becomes the root of the factor, an excluded term (cell line) is dropped.
- Subspecies are folded into their species: a label of three words or more becomes its first two words, with the id of that species in the datasets (e.g. *Mus musculus musculus* becomes *Mus musculus*, *Sus scrofa domesticus* becomes *Sus scrofa*).
- For assay, disease, sex and suspension type, `unknown`, `na` and an empty annotation are all "NA".
- The resources of a dataset are listed with a link to the dataset in each (several for datasets merged across resources), and datasets are flagged when merged across resources or annotated with several species.
- Datasets without a title (Single Cell Portal, ASAP, Bgee) get the collection name, or else the DOI label or the resource name, followed by the resource id.
- Each dataset gets one publication, so that the page can count publications instead of datasets. In this order: its DOI; else the DOI of the other datasets of its CELLxGENE collection; else its CELLxGENE collection; else the dataset itself (no publication information: the dataset counts as one study). Only CELLxGENE has real collections: in the scFAIR export, the collection id of Bgee and Single Cell Portal datasets is the dataset's own accession, and ASAP uses a placeholder (`ASAP000000`) for all projects. Result: 608 publications for 2,540 datasets (421 DOIs, 11 CELLxGENE collections without DOI, 176 datasets without publication information, mostly Single Cell Portal: only 8% of its datasets have a DOI in the export). The DOI also groups the datasets of one paper across resources (e.g. Li et al. 2022, the Fly Cell Atlas: 42 ASAP datasets and 1 Bgee dataset).

The terms file has one row per value of each field (`field`, `id`, `label`, `parent_ids`, `ancestor_ids`, `subsets`), including the publications (label: DOI label, e.g. "Siletti et al. (2023) Science", else collection name); for tissue, cell type and stage, the hierarchy comes from the TSV of step 4 and the subsets from the TSV of step 5. The page computes the counts and the roll-up in the browser.

## 7. Explore page

Static page (React, no server) with filters on the side and a pivot table (heatmap) of the number of datasets, e.g. species x tissue. Online at https://fbastian.github.io/BH2026-23-single-cell/, deployed on GitHub Pages from the `main` branch. See [explore/README.md](explore/README.md) to run it locally and for the deployment.
