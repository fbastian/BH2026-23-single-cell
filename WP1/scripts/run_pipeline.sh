#!/usr/bin/env bash
# Regenerate all WP1 data from scratch: steps 1 to 7 of WP1/README.md.
#
#   1. download the ontologies and the Bgee homology annotations into data/
#      (skipped with --skip-download)
#   2. merge the dataset exports            -> output/datasets_merged.json, merged_datasets.tsv
#   3. term remappings                      -> output/term_remappings.tsv
#   4. ancestor expansion                   -> output/term_ancestors.tsv
#   5. uninformative terms                  -> output/term_subsets.tsv
#   6. homology groups                      -> output/term_homology.tsv
#   7. export for the explore page          -> explore/ui/public/datasets.parquet, terms.parquet
#
# Inputs that are not downloaded (must be in data/): datasets_CxG.json, datasets_scFAIR.json,
# datasets_fenix.json. Hand-curated inputs: term_remappings_manual.tsv, term_subsets_manual.tsv.
#
# Usage (from anywhere; paths are relative to WP1/):
#   bash WP1/scripts/run_pipeline.sh [--skip-download]
# The Python interpreter is taken from $PYTHON (default: python, e.g. in the conda env bh23):
#   PYTHON=/path/to/envs/bh23/bin/python bash WP1/scripts/run_pipeline.sh
#
# Step 8 (the page) is not run: see WP1/explore/README.md (npm run build / npm run dev).

set -euo pipefail

SKIP_DOWNLOAD=0
for arg in "$@"; do
  case "$arg" in
    --skip-download) SKIP_DOWNLOAD=1 ;;
    -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg (see --help)" >&2; exit 2 ;;
  esac
done

PYTHON="${PYTHON:-python}"
cd "$(dirname "$0")/.."   # WP1/

step() { printf '\n=== %s\n' "$*"; }

need() {
  for f in "$@"; do
    [ -s "$f" ] || { echo "missing input file: WP1/$f" >&2; exit 1; }
  done
}

"$PYTHON" -c 'import pyarrow' 2>/dev/null \
  || { echo "$PYTHON cannot import pyarrow: use the bh23 environment (WP1/SETUP.md) or set PYTHON" >&2; exit 1; }

mkdir -p data output

step "1. Ontologies and homology annotations"
ONTOLOGIES=(uberon/composite-metazoan.obo uberon/collected-metazoan.obo cl.obo zfa.obo)
HOMOLOGY_URL=https://raw.githubusercontent.com/BgeeDB/anatomical-similarity-annotations/refs/heads/master/release/ancestral_taxa_homology_annotations.tsv
if [ "$SKIP_DOWNLOAD" -eq 1 ]; then
  echo "skipped (--skip-download), using the files in data/"
else
  for o in "${ONTOLOGIES[@]}"; do
    echo "downloading http://purl.obolibrary.org/obo/$o"
    curl -sSfL -o "data/$(basename "$o")" "http://purl.obolibrary.org/obo/$o"
  done
  echo "downloading $HOMOLOGY_URL"
  curl -sSfL -o "data/$(basename "$HOMOLOGY_URL")" "$HOMOLOGY_URL"
fi
for o in "${ONTOLOGIES[@]}"; do need "data/$(basename "$o")"; done
need "data/$(basename "$HOMOLOGY_URL")"
grep -m1 '^data-version:' data/composite-metazoan.obo || true

step "2. Merge the datasets"
need data/datasets_CxG.json data/datasets_scFAIR.json data/datasets_fenix.json
"$PYTHON" scripts/merge_datasets.py \
  --input cxg=data/datasets_CxG.json \
  --input scfair=data/datasets_scFAIR.json \
  --input fenix=data/datasets_fenix.json \
  --out output/datasets_merged.json \
  --report output/merged_datasets.tsv

step "3. Term remappings"
"$PYTHON" scripts/make_term_remappings.py \
  --composite data/composite-metazoan.obo \
  --datasets output/datasets_merged.json \
  --out output/term_remappings.tsv \
  --extra-obo data/collected-metazoan.obo data/cl.obo data/zfa.obo \
  --manual term_remappings_manual.tsv

step "4. Ancestor expansion"
"$PYTHON" scripts/make_term_ancestors.py \
  --composite data/composite-metazoan.obo \
  --datasets output/datasets_merged.json \
  --remappings output/term_remappings.tsv \
  --out output/term_ancestors.tsv

step "5. Uninformative terms"
"$PYTHON" scripts/make_term_subsets.py \
  --composite data/composite-metazoan.obo \
  --ancestors output/term_ancestors.tsv \
  --subsets non_informative upper_level ubprop:upper_level common_anatomy \
  --manual term_subsets_manual.tsv \
  --out output/term_subsets.tsv

step "6. Homology groups"
"$PYTHON" scripts/make_homology_groups.py \
  --annotations data/ancestral_taxa_homology_annotations.tsv \
  --composite data/composite-metazoan.obo \
  --ancestors output/term_ancestors.tsv \
  --out output/term_homology.tsv

step "7. Export for the explore page"
"$PYTHON" scripts/export_explore.py \
  --datasets output/datasets_merged.json \
  --remappings output/term_remappings.tsv \
  --ancestors output/term_ancestors.tsv \
  --subsets output/term_subsets.tsv \
  --homology output/term_homology.tsv \
  --out-datasets explore/ui/public/datasets.parquet \
  --out-terms explore/ui/public/terms.parquet

step "Done"
echo "Check the warnings above (unresolved remappings, terms outside their domain, ...)."
echo "Then rebuild or reload the page (WP1/explore/README.md) and commit the regenerated files."
