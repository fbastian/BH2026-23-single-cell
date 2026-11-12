import { useState } from "react";
import { Button, ListItemText, Menu, MenuItem } from "@mui/material";
import DownloadIcon from "@mui/icons-material/Download";
import {
  isPublication,
  matchingRows,
  termLabel,
  type ExperimentRow,
  type FacetField,
  type FacetSelection,
} from "../catalogQuery";

/** Ontology fields exported with both labels and ids. */
const TERM_FIELDS: FacetField[] = ["organism", "tissue", "cell_type", "development_stage", "assay", "disease", "sex"];
const SEP = "; ";
/** Column names of the exported fields (internal field ids otherwise). */
const COLUMN: Partial<Record<FacetField, string>> = { organism: "species", source: "resources" };
const col = (field: FacetField) => COLUMN[field] ?? field;

const clean = (v: string) => v.replace(/[\t\r\n]+/g, " ");
const labels = (field: FacetField, ids: string[]) => ids.map((id) => termLabel(field, id)).join(SEP);
const unique = (values: string[]) => [...new Set(values)];

function toTsv(header: string[], rows: string[][]): string {
  return [header, ...rows].map((r) => r.map((v) => clean(v ?? "")).join("\t")).join("\n") + "\n";
}

/** One row per dataset; multi-valued fields joined with "; ", term ids after remapping. */
export function datasetsTsv(rows: ExperimentRow[]): string {
  const header = [
    "dataset_id",
    "title",
    ...TERM_FIELDS.flatMap((f) => [col(f), `${col(f)}_ids`]),
    "suspension_type",
    "tissue_type",
    "cell_count",
    "resources",
    "publication",
    "publication_id",
    "publication_url",
    "flags",
    "explorer_url",
    "resource_urls",
    "h5ad_url",
  ];
  const body = rows.map((r) => {
    const pub = r.values.publication[0] ?? "";
    return [
      r.id,
      r.title,
      ...TERM_FIELDS.flatMap((f) => [labels(f, r.values[f]), r.values[f].join(SEP)]),
      labels("suspension_type", r.values.suspension_type),
      labels("tissue_type", r.values.tissue_type),
      String(r.cell_count),
      labels("source", r.values.source),
      isPublication(pub) ? termLabel("publication", pub) : "",
      pub,
      r.publication_url,
      labels("flags", r.values.flags),
      r.explorer_url,
      r.source_explorer_urls.join(SEP),
      r.h5ad_url,
    ];
  });
  return toTsv(header, body);
}

/**
 * One row per publication (or lone dataset without publication information) of the matching
 * datasets, with the union of their values, most datasets first.
 */
export function publicationsTsv(rows: ExperimentRow[]): string {
  const groups = new Map<string, ExperimentRow[]>();
  for (const r of rows) {
    const pub = r.values.publication[0] ?? r.id;
    const list = groups.get(pub);
    if (list) list.push(r);
    else groups.set(pub, [r]);
  }
  const fields: FacetField[] = ["organism", "tissue", "cell_type", "development_stage", "assay", "disease", "source"];
  const header = ["publication_id", "publication", "publication_url", "n_datasets", "total_cells", ...fields.map(col), "dataset_ids"];
  const body = [...groups.entries()]
    .sort((a, b) => b[1].length - a[1].length || a[0].localeCompare(b[0]))
    .map(([pub, ds]) => [
      pub,
      isPublication(pub) ? termLabel("publication", pub) : `${ds[0].title} (no publication information)`,
      ds[0].publication_url,
      String(ds.length),
      String(ds.reduce((sum, d) => sum + d.cell_count, 0)),
      ...fields.map((f) => labels(f, unique(ds.flatMap((d) => d.values[f])))),
      ds.map((d) => d.id).join(SEP),
    ]);
  return toTsv(header, body);
}

function download(filename: string, text: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: "text/tab-separated-values;charset=utf-8" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** "Export TSV" menu: the datasets matching the search and filters, or their publications. */
export function ExportButton({ query, facets }: { query: string; facets: FacetSelection }) {
  const [anchor, setAnchor] = useState<HTMLElement | null>(null);
  const date = new Date().toISOString().slice(0, 10);
  const run = async (kind: "datasets" | "publications") => {
    setAnchor(null);
    const rows = await matchingRows(query, facets);
    download(`${kind}_${date}.tsv`, kind === "datasets" ? datasetsTsv(rows) : publicationsTsv(rows));
  };
  return (
    <>
      <Button size="small" variant="outlined" startIcon={<DownloadIcon />} onClick={(e) => setAnchor(e.currentTarget)}>
        Export TSV
      </Button>
      <Menu anchorEl={anchor} open={Boolean(anchor)} onClose={() => setAnchor(null)}>
        <MenuItem onClick={() => run("datasets")}>
          <ListItemText primary="Datasets" secondary="one row per dataset, labels and ontology ids" />
        </MenuItem>
        <MenuItem onClick={() => run("publications")}>
          <ListItemText primary="Publications" secondary="one row per publication, with its datasets" />
        </MenuItem>
      </Menu>
    </>
  );
}
