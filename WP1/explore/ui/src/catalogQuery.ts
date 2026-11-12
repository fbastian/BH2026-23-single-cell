/**
 * In-browser catalog engine.
 *
 * Loads datasets.parquet and terms.parquet (written by WP1/scripts/export_explore.py) with
 * hyparquet and runs filter / facet / pivot in memory. Facet values are ids (ontology terms,
 * resource slugs, ...); labels come from the term table. Tissue, cell type and developmental
 * stage are hierarchical: a dataset annotated with a term also matches all its ancestors
 * (is_a + part_of within the factor domain, computed by make_term_ancestors.py).
 */
import { parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";

const BASE = import.meta.env.BASE_URL;
/** Cache-bust so GitHub Pages / browsers cannot keep a stale HTML 404 as parquet. */
const VERSION = "v=5";
export const DATASETS_URL = `${BASE}datasets.parquet?${VERSION}`;
export const TERMS_URL = `${BASE}terms.parquet?${VERSION}`;

export const HIERARCHICAL_FIELDS = ["tissue", "cell_type", "development_stage"] as const;

export const FACET_FIELDS = [
  "organism",
  "tissue",
  "cell_type",
  "development_stage",
  "assay",
  "disease",
  "sex",
  "suspension_type",
  "tissue_type",
  "source",
  "publication",
  "flags",
] as const;

export type FacetField = (typeof FACET_FIELDS)[number];
export type HierarchicalField = (typeof HIERARCHICAL_FIELDS)[number];
export type FacetSelection = Partial<Record<FacetField, string[]>>;

export const FIELD_TITLES: Record<FacetField, string> = {
  organism: "Species",
  tissue: "Tissue",
  cell_type: "Cell type",
  development_stage: "Developmental stage",
  assay: "Assay",
  disease: "Disease",
  sex: "Sex",
  suspension_type: "Suspension type",
  tissue_type: "Tissue type",
  source: "Resource",
  publication: "Publication",
  flags: "Flags",
};

/** Value id of 'not annotated' for the non-hierarchical fields (NA is the root for the others). */
export const NA = "na";

export type ExperimentRow = {
  objectID: string;
  id: string;
  title: string;
  cell_count: number;
  explorer_url: string;
  h5ad_url: string;
  doi: string;
  doi_label: string;
  collection_name: string;
  source_explorer_urls: string[];
  /** doi.org link or collection page of the publication; empty without publication information. */
  publication_url: string;
  search_blob: string;
  values: Record<FacetField, string[]>;
  /** Annotated terms and all their ancestors, per hierarchical field. */
  expanded: Record<HierarchicalField, Set<string>>;
};

export type Term = {
  id: string;
  label: string;
  parents: string[];
  ancestors: string[];
  /** Subsets of the term (e.g. upper_level): hidden in filter trees and roll-up columns. */
  subsets: string[];
};

type Catalog = {
  rows: ExperimentRow[];
  terms: Record<FacetField, Map<string, Term>>;
  children: Record<HierarchicalField, Map<string, string[]>>;
  roots: Record<HierarchicalField, string[]>;
  /** Datasets per publication id. */
  publicationSizes: Map<string, number>;
};

export type FacetCount = { value: string; count: number };

/**
 * One pivot axis: a field; for hierarchical fields a roll-up node (null = raw terms); and the
 * filter selection on that field (scope), so that the axis only shows the selected values
 * (and, for hierarchical fields, their descendants).
 */
export type PivotAxis = { field: FacetField; node: string | null; scope: string[] };

/** Roll-up node whose children are the selected filter terms (several terms selected). */
export const SELECTED_NODE = "__selected__";

export function stripExact(value: string): string {
  return value.startsWith(EXACT_PREFIX) ? value.slice(EXACT_PREFIX.length) : value;
}

/** Pivot value for datasets annotated with the roll-up node itself, not with a child. */
export const SELF_PREFIX = "__self__:";

/**
 * Filter value matching datasets annotated exactly with a term, not with its descendants
 * (used when clicking a raw pivot cell, so the listed datasets are those of the cell).
 */
export const EXACT_PREFIX = "=";

export type PivotResult = {
  rows: string[];
  cols: string[];
  /** Counted units (datasets or publications) per cell. */
  matrix: number[][];
  /** Datasets per cell (same as matrix when counting datasets). */
  datasetMatrix: number[][];
  rowTotals: number[];
  colTotals: number[];
  maxCell: number;
  /** Datasets in the pivot. */
  found: number;
  /** Counted units in the pivot (datasets or publications). */
  units: number;
};

/** Publication ids of datasets without publication information (one dataset = one study). */
export const DATASET_PUBLICATION_PREFIX = "dataset:";

/** Whether a publication id stands for a real publication or collection, not a lone dataset. */
export function isPublication(id: string | undefined): id is string {
  return Boolean(id) && !id!.startsWith(DATASET_PUBLICATION_PREFIX);
}

/** Number of datasets of a publication (whole catalog). */
export function publicationSize(id: string): number {
  return catalog?.publicationSizes.get(id) ?? 0;
}

let catalog: Catalog | null = null;
let catalogPromise: Promise<Catalog> | null = null;

export function isHierarchical(field: string): field is HierarchicalField {
  return (HIERARCHICAL_FIELDS as readonly string[]).includes(field);
}

function asStrings(value: unknown): string[] {
  if (value == null) return [];
  if (Array.isArray(value)) return value.map(String).filter((v) => v.length > 0);
  return [String(value)];
}

function assertParquetMagic(buffer: ArrayBuffer, url: string): void {
  const bytes = new Uint8Array(buffer);
  const text = (from: number) => String.fromCharCode(...bytes.slice(from, from + 4));
  if (bytes.length < 8 || text(0) !== "PAR1" || text(bytes.length - 4) !== "PAR1") {
    const sniff = new TextDecoder().decode(bytes.slice(0, 48)).replace(/\s+/g, " ");
    throw new Error(
      `Invalid parquet from ${url} (starts ${JSON.stringify(sniff)}). Hard-refresh the page, or re-export: python scripts/export_explore.py`
    );
  }
}

async function readParquet(url: string): Promise<Record<string, unknown>[]> {
  // Full-file fetch (not range/URL mode): GitHub Pages gzips .parquet and
  // applies byte ranges to the gzip stream, which breaks hyparquet footers.
  // no-store avoids 304 reuse of a previously cached HTML/error body.
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(
      `Failed to fetch ${url} (${response.status}). Run: python scripts/export_explore.py`
    );
  }
  const buffer = await response.arrayBuffer();
  assertParquetMagic(buffer, url);
  return (await parquetReadObjects({ file: buffer, compressors })) as Record<string, unknown>[];
}

function buildCatalog(
  datasetRows: Record<string, unknown>[],
  termRows: Record<string, unknown>[]
): Catalog {
  const terms = Object.fromEntries(
    FACET_FIELDS.map((f) => [f, new Map<string, Term>()])
  ) as Catalog["terms"];
  for (const r of termRows) {
    const field = String(r.field) as FacetField;
    terms[field]?.set(String(r.id), {
      id: String(r.id),
      label: String(r.label),
      parents: asStrings(r.parent_ids),
      ancestors: asStrings(r.ancestor_ids),
      subsets: asStrings(r.subsets),
    });
  }

  const children = {} as Catalog["children"];
  const roots = {} as Catalog["roots"];
  for (const field of HIERARCHICAL_FIELDS) {
    const kids = new Map<string, string[]>();
    const top: string[] = [];
    for (const term of terms[field].values()) {
      if (!term.parents.length) top.push(term.id);
      for (const p of term.parents) {
        const list = kids.get(p);
        if (list) list.push(term.id);
        else kids.set(p, [term.id]);
      }
    }
    children[field] = kids;
    roots[field] = top;
  }

  const rows = datasetRows.map((r): ExperimentRow => {
    const id = String(r.id ?? "");
    const values = Object.fromEntries(
      FACET_FIELDS.map((f) => [f, asStrings(r[f])])
    ) as Record<FacetField, string[]>;
    const expanded = Object.fromEntries(
      HIERARCHICAL_FIELDS.map((f) => {
        const set = new Set<string>();
        for (const t of values[f]) {
          set.add(t);
          for (const a of terms[f].get(t)?.ancestors ?? []) set.add(a);
        }
        return [f, set];
      })
    ) as Record<HierarchicalField, Set<string>>;
    return {
      objectID: id,
      id,
      title: String(r.title ?? id),
      cell_count: Number(r.cell_count ?? 0),
      explorer_url: String(r.explorer_url ?? ""),
      h5ad_url: String(r.h5ad_url ?? ""),
      doi: String(r.doi ?? ""),
      doi_label: String(r.doi_label ?? ""),
      collection_name: String(r.collection_name ?? ""),
      source_explorer_urls: asStrings(r.source_explorer_urls),
      publication_url: String(r.publication_url ?? ""),
      search_blob: String(r.search_blob ?? ""),
      values,
      expanded,
    };
  });
  const publicationSizes = new Map<string, number>();
  for (const row of rows) {
    for (const p of row.values.publication) publicationSizes.set(p, (publicationSizes.get(p) ?? 0) + 1);
  }
  return { rows, terms, children, roots, publicationSizes };
}

async function loadCatalog(): Promise<Catalog> {
  if (!catalogPromise) {
    catalogPromise = Promise.all([readParquet(DATASETS_URL), readParquet(TERMS_URL)])
      .then(([datasetRows, termRows]) => {
        catalog = buildCatalog(datasetRows, termRows);
        return catalog;
      })
      .catch((err) => {
        catalogPromise = null;
        throw err;
      });
  }
  return catalogPromise;
}

/** Loads the catalog (both Parquet files); resolves when the page can render. */
export async function initCatalog(): Promise<void> {
  await loadCatalog();
}

function loaded(): Catalog {
  if (!catalog) throw new Error("catalog not loaded");
  return catalog;
}

/** Label of a value id (synchronous: the catalog is loaded before anything renders). */
export function termLabel(field: FacetField, id: string): string {
  if (id === SELECTED_NODE) return "selected terms";
  if (id.startsWith(SELF_PREFIX)) {
    return `${termLabel(field, id.slice(SELF_PREFIX.length))} (not more specific)`;
  }
  if (id.startsWith(EXACT_PREFIX)) {
    return `${termLabel(field, id.slice(EXACT_PREFIX.length))} (exactly)`;
  }
  return catalog?.terms[field].get(id)?.label ?? id;
}

export function labelsOf(field: FacetField, ids: string[]): string[] {
  return ids.map((id) => termLabel(field, id));
}

/**
 * Whether a term is hidden in the filter trees and roll-up columns: it is in one of the
 * subsets of the term subset TSV (e.g. upper_level, non_informative).
 */
export function isHidden(field: FacetField, id: string): boolean {
  return Boolean(catalog?.terms[field].get(id)?.subsets.length);
}

const visibleChildrenCache = new Map<string, string[]>();

/** Children of a term, hidden children being replaced by their own (visible) children. */
export function childrenOf(field: HierarchicalField, id: string): string[] {
  const key = `${field}\0${id}`;
  let out = visibleChildrenCache.get(key);
  if (!out) {
    const found = new Set<string>();
    for (const c of loaded().children[field].get(id) ?? []) {
      if (isHidden(field, c)) childrenOf(field, c).forEach((g) => found.add(g));
      else found.add(c);
    }
    out = [...found];
    visibleChildrenCache.set(key, out);
  }
  return out;
}

export function rootsOf(field: HierarchicalField): string[] {
  return loaded().roots[field];
}

/** Top terms of a filter tree: the roots, or their visible children for hidden roots. */
export function topTermsOf(field: HierarchicalField): string[] {
  return [...new Set(rootsOf(field).flatMap((r) => (isHidden(field, r) ? childrenOf(field, r) : [r])))];
}

/**
 * A path from a root down to the term (following the first parent at each step), without
 * the hidden terms (the root is kept as a starting point).
 */
export function pathTo(field: HierarchicalField, id: string): string[] {
  const terms = loaded().terms[field];
  const path = [id];
  for (let t = terms.get(id); t?.parents.length && path.length < 100; t = terms.get(t.parents[0])) {
    path.unshift(t.parents[0]);
  }
  return path.filter((t, i) => i === 0 || i === path.length - 1 || !isHidden(field, t));
}

/** All value ids of a field (for hierarchical fields: annotated terms and their ancestors). */
export function allValues(field: FacetField): string[] {
  return [...loaded().terms[field].keys()];
}

/** Value ids of a field that are not hidden. */
export function visibleValues(field: FacetField): string[] {
  return allValues(field).filter((id) => !isHidden(field, id));
}

/** Values a filter on the field tests: with ancestors for the hierarchical fields. */
function matchValues(row: ExperimentRow, field: FacetField): Set<string> | string[] {
  return isHierarchical(field) ? row.expanded[field] : row.values[field];
}

function has(values: Set<string> | string[], v: string): boolean {
  return values instanceof Set ? values.has(v) : values.includes(v);
}

function rowMatchesFacets(row: ExperimentRow, facets: FacetSelection, omitFacet?: string) {
  for (const [field, wanted] of Object.entries(facets) as [FacetField, string[]][]) {
    if (!wanted?.length || field === omitFacet) continue;
    const present = matchValues(row, field);
    const hit = (v: string) =>
      v.startsWith(EXACT_PREFIX)
        ? row.values[field].includes(v.slice(EXACT_PREFIX.length))
        : has(present, v);
    if (!wanted.some(hit)) return false;
  }
  return true;
}

function rowMatchesQuery(row: ExperimentRow, query: string): boolean {
  const q = query.trim().toLowerCase();
  if (!q) return true;
  return row.search_blob.toLowerCase().includes(q) || row.id.toLowerCase().includes(q);
}

function filteredRows(rows: ExperimentRow[], query: string, facets: FacetSelection, omitFacet?: string) {
  return rows.filter((row) => rowMatchesQuery(row, query) && rowMatchesFacets(row, facets, omitFacet));
}

/** Publication (or lone dataset) of a dataset: the unit counted in "publications" mode. */
function publicationOf(row: ExperimentRow): string {
  return row.values.publication[0] ?? row.id;
}

/**
 * Number of datasets (or of distinct publications) per value of a field, with the filters of
 * the other fields. For the hierarchical fields, a dataset counts for its terms and all their
 * ancestors.
 */
export async function facetCounts(
  field: FacetField,
  query: string,
  facets: FacetSelection,
  countBy: CountUnit = "datasets"
): Promise<Map<string, number>> {
  const rows = filteredRows((await loadCatalog()).rows, query, facets, field);
  const units = new Map<string, Set<string>>();
  for (const row of rows) {
    const unit = countBy === "publications" ? publicationOf(row) : row.id;
    for (const value of new Set(matchValues(row, field))) {
      const set = units.get(value);
      if (set) set.add(unit);
      else units.set(value, new Set([unit]));
    }
  }
  return new Map([...units].map(([value, set]) => [value, set.size]));
}

export async function searchRows(
  query: string,
  facets: FacetSelection,
  options: { offset: number; limit: number }
): Promise<{ rows: ExperimentRow[]; found: number; publications: number }> {
  const { offset, limit } = options;
  const matched = filteredRows((await loadCatalog()).rows, query, facets).sort(
    (a, b) => b.cell_count - a.cell_count || a.title.localeCompare(b.title)
  );
  return {
    found: matched.length,
    publications: new Set(matched.map(publicationOf)).size,
    rows: matched.slice(offset, offset + limit),
  };
}

/** Whether a value is in the filter selection of its field (or under a selected term). */
function inScope(field: FacetField, value: string, scope: string[]): boolean {
  if (!scope.length) return true;
  return scope.some((s) => {
    if (s.startsWith(EXACT_PREFIX)) return value === stripExact(s);
    if (value === s) return true;
    return isHierarchical(field) && (loaded().terms[field].get(value)?.ancestors.includes(s) ?? false);
  });
}

/**
 * Values of a dataset on a pivot axis. Raw: its (remapped) annotations within the filter
 * selection of the field. Roll-up on a node: all terms under the node the dataset counts
 * for (see rollupValues). Roll-up on SELECTED_NODE: the roll-up under each selected term,
 * merged.
 */
function axisValues(row: ExperimentRow, axis: PivotAxis): string[] {
  if (!isHierarchical(axis.field) || axis.node == null) {
    return row.values[axis.field].filter((v) => inScope(axis.field, v, axis.scope));
  }
  const expanded = row.expanded[axis.field];
  if (axis.node === SELECTED_NODE) {
    // Several selected terms: the roll-up under each of them, merged.
    return axis.scope.flatMap((s) => {
      if (s.startsWith(EXACT_PREFIX)) {
        const id = stripExact(s);
        return row.values[axis.field].includes(id) ? [SELF_PREFIX + id] : [];
      }
      return rollupValues(expanded, axis.field as HierarchicalField, s);
    });
  }
  return rollupValues(expanded, axis.field, axis.node);
}

/**
 * All terms under the node (at any depth) that the dataset is annotated with or under, or
 * SELF_PREFIX+node if it is annotated with the node itself and nothing more specific, or
 * nothing if it is outside the node.
 */
function rollupValues(expanded: Set<string>, field: HierarchicalField, node: string): string[] {
  if (!expanded.has(node)) return [];
  const below = [...expanded].filter((t) => !isHidden(field, t) && ancestorSet(field, t).has(node));
  return below.length ? below : [SELF_PREFIX + node];
}

const ancestorSets = new Map<string, Set<string>>();

/** Ancestors of a term as a set (cached). */
function ancestorSet(field: HierarchicalField, id: string): Set<string> {
  const key = `${field}\0${id}`;
  let set = ancestorSets.get(key);
  if (!set) {
    set = new Set(loaded().terms[field].get(id)?.ancestors ?? []);
    ancestorSets.set(key, set);
  }
  return set;
}

/** Pivot counts: datasets, or publications (each dataset belongs to one publication or study). */
export type CountUnit = "datasets" | "publications";

export async function pivotCounts(
  rowAxis: PivotAxis,
  colAxis: PivotAxis,
  query: string,
  facets: FacetSelection,
  options: { maxRows?: number; maxCols?: number; countBy?: CountUnit } = {}
): Promise<PivotResult> {
  const maxRows = options.maxRows ?? 40;
  const maxCols = options.maxCols ?? 30;
  const byPublication = options.countBy === "publications";
  const matched = filteredRows((await loadCatalog()).rows, query, facets);
  // Distinct units (dataset ids, or publication ids) per cell, row value and column value.
  const pairUnits = new Map<string, Set<string>>();
  const pairDatasets = new Map<string, number>();
  const rowUnits = new Map<string, Set<string>>();
  const colUnits = new Map<string, Set<string>>();
  const allUnits = new Set<string>();
  const add = (map: Map<string, Set<string>>, key: string, unit: string) => {
    const set = map.get(key);
    if (set) set.add(unit);
    else map.set(key, new Set([unit]));
  };
  let found = 0;

  for (const row of matched) {
    const rows = [...new Set(axisValues(row, rowAxis))];
    const cols = [...new Set(axisValues(row, colAxis))];
    if (!rows.length || !cols.length) continue;
    found += 1;
    const unit = byPublication ? publicationOf(row) : row.id;
    allUnits.add(unit);
    // A dataset (or publication) counts once per row value, column value and cell.
    for (const r of rows) add(rowUnits, r, unit);
    for (const c of cols) add(colUnits, c, unit);
    for (const r of rows) {
      for (const c of cols) {
        const key = `${r}\0${c}`;
        add(pairUnits, key, unit);
        pairDatasets.set(key, (pairDatasets.get(key) ?? 0) + 1);
      }
    }
  }

  const top = (units: Map<string, Set<string>>, max: number) =>
    [...units.entries()]
      .sort((a, b) => b[1].size - a[1].size || a[0].localeCompare(b[0]))
      .slice(0, max)
      .map(([value]) => value);
  const rowOrder = top(rowUnits, maxRows);
  const colOrder = top(colUnits, maxCols);

  let maxCell = 0;
  const matrix = rowOrder.map((r) =>
    colOrder.map((c) => {
      const n = pairUnits.get(`${r}\0${c}`)?.size ?? 0;
      if (n > maxCell) maxCell = n;
      return n;
    })
  );

  return {
    rows: rowOrder,
    cols: colOrder,
    matrix,
    datasetMatrix: rowOrder.map((r) => colOrder.map((c) => pairDatasets.get(`${r}\0${c}`) ?? 0)),
    rowTotals: rowOrder.map((r) => rowUnits.get(r)?.size ?? 0),
    colTotals: colOrder.map((c) => colUnits.get(c)?.size ?? 0),
    maxCell,
    found,
    units: allUnits.size,
  };
}
