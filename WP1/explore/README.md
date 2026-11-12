# Explore page: species x factor table

**Online: https://fbastian.github.io/BH2026-23-single-cell/**

Static page (React + Vite + MUI, Nivo heatmap, [hyparquet](https://github.com/hyparam/hyparquet)) to browse the merged datasets of the WP1 pipeline. It needs no server: the browser loads two Parquet files and computes the filters, counts and ontology roll-up in memory.

| Piece | Role |
|-------|------|
| [`../scripts/export_explore.py`](../scripts/export_explore.py) | step 7 of the [WP1 pipeline](../README.md): writes `ui/public/datasets.parquet` (one row per dataset) and `ui/public/terms.parquet` (labels and hierarchy of the values) |
| `ui/src/catalogQuery.ts` | loads the Parquet files; filters, facet counts, pivot counts, roll-up |
| `ui/src/components/` | filters (`FacetList`), pivot heatmap (`PivotHeatmap`), dataset cards and table |

## What the page shows

- **Filters** (left): species, tissue, cell type, developmental stage, assay, disease, sex, suspension type, tissue type, resource, publication (searchable), flags (merged across resources, several species). Values of a field are combined with OR, fields with AND. Tissue, cell type and stage are trees (composite-metazoan Uberon and CL, `is_a` + `part_of` within the factor domain): the arrow shows the more specific terms, and ticking a term also matches the datasets annotated with any of its descendants. The search box of a tree searches all its terms.
- **Pivot** (default view, opening on species x tissue in roll-up mode): number of datasets for two fields, e.g. species x tissue. For tissue, cell type and stage, each axis is either:
  - **Raw**: the terms as annotated (after remapping); a dataset annotated with "hippocampal formation" counts only in that column;
  - **Roll-up**: a dataset also counts for all the more general terms of its annotations (the hippocampal formation dataset also counts for brain, forebrain, ...). The columns are all the terms under the current term (the root at first), at any depth, sorted by number of datasets ("Show top" keeps the first ones). Click an underlined label to show only the terms under it, use the breadcrumb to go back up, or jump to any term. The column named after the current term itself holds the datasets annotated with that term and with no term below it (they would otherwise appear in no column; its tooltip says so): at the root, these are the datasets without annotation.

  Terms of the subsets tagged in step 5 of the pipeline (Uberon `non_informative`, `upper_level` and `common_anatomy`, and hand-curated terms; e.g. regional part of brain, organ, material anatomical entity) are hidden in the filter trees and the roll-up columns: their more specific terms are shown instead (they still count, through these). In raw mode, terms annotated to datasets are shown even if they are in these subsets.

  The axes follow the side filters of their field: with "brain" ticked in the tissue filter, the tissue roll-up starts under brain (the breadcrumb still goes up); with several tissues ticked, the columns are the terms under each of them, merged (e.g. the parts of brain and of heart); in raw mode, and for the fields without hierarchy (e.g. species), only the selected values (and their more specific terms) are shown, not the other annotations of the matching datasets.

  **Homologous terms merged into one column** (columns only, in raw and roll-up modes). When terms of the column field are ticked in the side filters, the pivot also includes their homologous terms, and a notice names them: ticking lung brings swim bladder (and the zebrafish datasets), so the merged column appears without knowing the homology beforehand, as soon as both have data with the other filters (with only human ticked, the column stays "lung"). This applies to the pivot only: the side filters, cards, table and export keep the literal selection. When at least two terms of a homology group (step 6 of the pipeline, from the Bgee anatomical similarity annotations, e.g. lung and swim bladder, homologous in Gnathostomata) are columns, they form a single column, labelled with the terms of the group joined by `|` in italics (e.g. under "pectoral appendage": "pectoral fin | arm | forelimb"); the columns named after ticked terms themselves are merged too (e.g. with lung and swim bladder ticked, the column "lung | swim bladder" holds the datasets annotated with lung or swim bladder and with no term below them); it counts each dataset annotated with any of them once, and its tooltip gives the taxa in which the homology holds. Its "list" button lists the datasets of all its terms. Rows, filters, cards, table and export are not affected.

  **Count: Datasets / Publications** (switch in the top bar, next to the views; it applies to the side filters, the pivot and the result summary). With "Publications", each number counts the distinct publications of the datasets: the datasets of a publication (same DOI, or same CELLxGENE collection without DOI) count once, and a dataset without publication information counts as one study (see step 7 of the [pipeline](../README.md)). The tooltip of a cell also gives its number of datasets.

  A dataset counts in every row and column it is annotated with: a dataset with two species counts in both species rows, a dataset with five tissues in five columns. Click a cell to list its datasets, or the "list" button under a column to list all the datasets of the column. This sets the filter of the field to that value (exact for raw terms), replacing the field's current selection: values of a field combine with OR, so keeping "brain" next to "forebrain" would list all brain datasets.
- **Cards / Table**: the datasets matching the filters, with all the fields of the side filters (species, tissue, cell type, developmental stage, assay, disease, sex, suspension type, tissue type, resource, publication, flags), links to each resource, to the publication (DOI or CELLxGENE collection; its tooltip gives its number of datasets) and to the H5AD file. Values are shown as chips whose text wraps (never cut); long lists show their first values and a "+N more" chip that unfolds the whole list. The table has the title, publication, species and number of cells first, then one column per other field (hide some with "Columns"); it scrolls horizontally with the title column pinned, and can be grouped by any field, e.g. publication.

- **Export TSV** (cards and table views): all the datasets matching the search and filters (not only the displayed page), as a tab-separated file built in the browser:
  - *Datasets*: one row per dataset: `dataset_id`, `title`, `publication`, `publication_id`, `publication_url`, then label and ids of species, tissue, cell type, developmental stage, assay, disease, sex (e.g. `tissue`, `tissue_ids`; terms after remapping to composite-metazoan Uberon/CL), `suspension_type`, `tissue_type`, `cell_count`, `resources`, `flags`, `explorer_url`, `resource_urls`, `h5ad_url`. Multiple values are separated by `; `.
  - *Publications*: one row per publication of these datasets (or per dataset without publication information), most datasets first: `publication_id`, `publication`, `publication_url`, `n_datasets`, `total_cells`, the values of species, tissue, cell type, developmental stage, assay, disease and resources over its datasets, and `dataset_ids`.

**Shareable links and back button.** The address of the page always reflects its state, so it can be copied and shared (or use the link button in the top bar); opening it restores the same view. Each change (view, filters, count, page, pivot axes or drill-down, ...) is a step of the browser history, so the back and forward buttons go through the previous states (typing in the search box is not a step; "Clear" is, so back undoes it). Parameters: `view` (pivot, cards, table), `count` (datasets, publications), `q` (search), one parameter per filtered field with comma-separated value ids of the terms file (`organism`, `tissue`, `cell_type`, `development_stage`, `assay`, `disease`, `sex`, `suspension_type`, `tissue_type`, `source`, `publication`, `flags`; `=` before an id means exactly this term), `page`, `group` (table grouping), and for the pivot `rows`, `cols` (fields), `rows_mode`, `cols_mode` (raw, rollup), `rows_path`, `cols_path` (drill-down position) and `top`. Example: `?view=pivot&count=publications&tissue=UBERON:0002048&organism=NCBITaxon:9606,NCBITaxon:7955&cols=tissue&cols_mode=rollup`. Personal display choices (table columns, open filter panels) stay in the browser and are not part of the link.

## Run locally

Prerequisite: Node.js 20 or later. From `WP1/`, after steps 1 to 7 of the [pipeline](../README.md):

```bash
cd explore/ui
npm ci
npm run dev        # http://localhost:5173, reloads on code changes
# or the production build:
npm run build && npm run preview
```

After re-running the export, hard-refresh the browser.

## GitHub Pages

[`.github/workflows/deploy-explore-pages.yml`](../../.github/workflows/deploy-explore-pages.yml) builds `ui/` and pushes it to the `gh-pages` branch on every push to `main` that changes `WP1/explore/ui/` (or by hand: Actions, "Deploy explore UI to GitHub Pages", Run workflow). It does not run the Python export: commit the regenerated Parquet files in `ui/public/`.

GitHub Pages is enabled for this repository: the page is at https://fbastian.github.io/BH2026-23-single-cell/. After a push, the new version is online once the workflow and the "pages build and deployment" run have finished (Actions tab, about a minute); hard-refresh the browser.

In a fork, one-time setup: Settings, Pages, Build and deployment, Source "Deploy from a branch", branch `gh-pages`, folder `/ (root)`. The page is then at `https://<owner>.github.io/BH2026-23-single-cell/` (the build uses `VITE_BASE_PATH=/BH2026-23-single-cell/`).
