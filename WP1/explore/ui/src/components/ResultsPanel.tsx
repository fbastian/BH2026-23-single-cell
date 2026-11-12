import { useEffect, useMemo, useState, type MouseEvent } from "react";
import {
  Button,
  Checkbox,
  FormControl,
  InputLabel,
  Link,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Select,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Tooltip,
  Typography,
} from "@mui/material";
import ViewColumnIcon from "@mui/icons-material/ViewColumn";
import {
  EXACT_PREFIX,
  FACET_FIELDS,
  FIELD_TITLES,
  isHierarchical,
  labelsOf,
  type ExperimentRow,
  type FacetField,
  type FacetSelection,
} from "../catalogQuery";
import { EllipsisChip } from "./EllipsisChip";
import { GROUP_BY_OPTIONS, groupHits, type GroupByField } from "./groupBy";
import {
  FlagChips,
  Hit,
  PublicationChip,
  SourceChips,
  ValueChips,
  formatCompactCount,
  hitTitle,
  publicationLabel,
  type ExperimentHit,
} from "./Hit";
import { PivotHeatmap } from "./PivotHeatmap";
import type { ResultsView } from "./ViewToggle";

const GROUP_STORAGE_KEY = "wp1-explore-group-by";
// v3: all the side-filter fields as columns, shown by default (2026-11-12)
const COLUMNS_STORAGE_KEY = "wp1-explore-table-columns-v3";

type TableColumnId = "title" | "cells" | FacetField;

/** Column widths, wide enough for whole words in the chips (e.g. CELLxGENE, annotated). */
const COLUMN_WIDTHS: Record<FacetField, number> = {
  organism: 130,
  tissue: 190,
  cell_type: 210,
  development_stage: 170,
  assay: 130,
  disease: 150,
  sex: 110,
  suspension_type: 115,
  tissue_type: 125,
  source: 125,
  publication: 190,
  flags: 125,
};

const fieldColumn = (field: FacetField) => ({
  id: field as TableColumnId,
  label: FIELD_TITLES[field],
  width: COLUMN_WIDTHS[field],
});

/** Columns: the title, species, the number of cells, then the other side-filter fields. */
const TABLE_COLUMNS: { id: TableColumnId; label: string; width: number; required?: boolean }[] = [
  { id: "title", label: "Title", width: 260, required: true },
  fieldColumn("organism"),
  { id: "cells", label: "Cells", width: 72 },
  ...FACET_FIELDS.filter((field) => field !== "organism").map(fieldColumn),
];

const DEFAULT_COLUMNS: TableColumnId[] = TABLE_COLUMNS.map((col) => col.id);

/** Opaque table background: the theme's paper colour is translucent, which would let the
 * columns scrolled under the pinned title column show through. */
const TABLE_BG = "#ffffff";

/** The title column stays visible when the table scrolls horizontally. */
const STICKY_SX = {
  position: "sticky",
  left: 0,
  zIndex: 1,
  bgcolor: TABLE_BG,
  boxShadow: "inset -1px 0 0 rgba(0, 0, 0, 0.08)",
} as const;

function loadStoredGroupBy(): GroupByField {
  try {
    const value = localStorage.getItem(GROUP_STORAGE_KEY);
    if (GROUP_BY_OPTIONS.some((option) => option.value === value)) return value as GroupByField;
  } catch {
    // ignore
  }
  return "none";
}

function loadStoredColumns(): TableColumnId[] {
  try {
    const parsed = JSON.parse(localStorage.getItem(COLUMNS_STORAGE_KEY) ?? "null") as unknown;
    if (!Array.isArray(parsed)) return DEFAULT_COLUMNS;
    const valid = parsed.filter((id): id is TableColumnId => TABLE_COLUMNS.some((col) => col.id === id));
    if (!valid.includes("title")) valid.unshift("title");
    return valid;
  } catch {
    return DEFAULT_COLUMNS;
  }
}

function groupDomId(key: string): string {
  const slug = key.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 80);
  return `group-section-${slug || "none"}`;
}

type Props = {
  view: ResultsView;
  rows: ExperimentRow[];
  query: string;
  facets: FacetSelection;
  onOpenHits: () => void;
  onApplyFacet: (field: FacetField, value: string) => void;
};

function HitRow({ hit, columns }: { hit: ExperimentHit; columns: typeof TABLE_COLUMNS }) {
  const title = hitTitle(hit);
  return (
    <TableRow hover>
      {columns.map((col) => {
        const sx = { verticalAlign: "top", px: 1, ...(col.id === "title" ? STICKY_SX : {}) };
        if (col.id === "title") {
          return (
            <TableCell key={col.id} sx={sx}>
              {hit.explorer_url ? (
                <Link href={hit.explorer_url} target="_blank" rel="noreferrer" underline="hover">
                  {title}
                </Link>
              ) : (
                title
              )}
            </TableCell>
          );
        }
        if (col.id === "source") {
          return (
            <TableCell key={col.id} sx={sx}>
              <Stack direction="row" flexWrap="wrap" useFlexGap gap={0.5}>
                <SourceChips hit={hit} />
              </Stack>
            </TableCell>
          );
        }
        if (col.id === "publication" || col.id === "flags") {
          const empty = col.id === "publication" ? !publicationLabel(hit) : !hit.values.flags.length;
          return (
            <TableCell key={col.id} sx={sx}>
              {empty ? (
                <Typography variant="body2" color="text.secondary">
                  —
                </Typography>
              ) : (
                <Stack direction="row" flexWrap="wrap" useFlexGap gap={0.5}>
                  {col.id === "publication" ? <PublicationChip hit={hit} /> : <FlagChips hit={hit} />}
                </Stack>
              )}
            </TableCell>
          );
        }
        if (col.id === "cells") {
          return (
            <TableCell key={col.id} align="right" sx={sx}>
              <Tooltip title={hit.cell_count.toLocaleString()} enterDelay={400}>
                <span>{formatCompactCount(hit.cell_count)}</span>
              </Tooltip>
            </TableCell>
          );
        }
        const field = col.id as FacetField;
        return (
          <TableCell key={col.id} sx={sx}>
            <ValueChips values={labelsOf(field, hit.values[field])} limit={4} />
          </TableCell>
        );
      })}
    </TableRow>
  );
}

export function ResultsPanel({ view, rows, query, facets, onOpenHits, onApplyFacet }: Props) {
  const [groupBy, setGroupBy] = useState<GroupByField>(() => loadStoredGroupBy());
  const [columns, setColumns] = useState<TableColumnId[]>(() => loadStoredColumns());
  const [activeChip, setActiveChip] = useState<string | null>(null);
  const [columnsMenuAnchor, setColumnsMenuAnchor] = useState<HTMLElement | null>(null);

  useEffect(() => {
    try {
      localStorage.setItem(GROUP_STORAGE_KEY, groupBy);
    } catch {
      // ignore
    }
    setActiveChip(null);
  }, [groupBy]);

  useEffect(() => {
    try {
      localStorage.setItem(COLUMNS_STORAGE_KEY, JSON.stringify(columns));
    } catch {
      // ignore
    }
  }, [columns]);

  const visibleColumns = useMemo(() => TABLE_COLUMNS.filter((col) => columns.includes(col.id)), [columns]);
  const tableWidth = visibleColumns.reduce((sum, col) => sum + col.width, 0);

  const groups = useMemo(() => (view === "table" ? groupHits(rows, groupBy) : []), [view, rows, groupBy]);

  const topGroups = useMemo(() => {
    if (groupBy === "none") return [];
    return [...groups]
      .sort((a, b) => b.hits.length - a.hits.length || a.label.localeCompare(b.label, undefined, { sensitivity: "base" }))
      .slice(0, 10);
  }, [groups, groupBy]);

  const focusGroup = (key: string) => {
    setActiveChip(key);
    if (groupBy === "none") return;
    // Groups are the annotated values: filter on exactly that term, not its descendants.
    onApplyFacet(groupBy, isHierarchical(groupBy) ? EXACT_PREFIX + key : key);
    document.getElementById(groupDomId(key))?.scrollIntoView({ behavior: "smooth", block: "start" });
  };

  const toggleColumn = (id: TableColumnId) => {
    if (id === "title") return;
    setColumns((prev) => (prev.includes(id) ? prev.filter((col) => col !== id) : [...prev, id]));
  };

  if (view === "pivot") {
    return (
      <PivotHeatmap
        query={query}
        facets={facets}
        onApplyCell={(cells) => {
          for (const { field, value } of cells) onApplyFacet(field, value);
          onOpenHits();
        }}
      />
    );
  }

  if (view === "cards") {
    return (
      <Stack
        direction="row"
        flexWrap="wrap"
        useFlexGap
        gap={1.25}
        sx={{
          "& > *": {
            flex: { xs: "1 1 100%", sm: "1 1 calc(50% - 10px)", lg: "1 1 calc(33.333% - 12px)" },
            maxWidth: { sm: "calc(50% - 10px)", lg: "calc(33.333% - 12px)" },
          },
        }}
      >
        {rows.map((hit) => (
          <Hit key={hit.objectID} hit={hit} />
        ))}
      </Stack>
    );
  }

  return (
    <Stack spacing={2}>
      <Stack direction={{ xs: "column", sm: "row" }} spacing={1.5} alignItems={{ sm: "center" }} flexWrap="wrap" useFlexGap>
        <FormControl size="small" sx={{ minWidth: 180 }}>
          <InputLabel id="group-by-label">Group by</InputLabel>
          <Select
            labelId="group-by-label"
            label="Group by"
            value={groupBy}
            onChange={(event) => setGroupBy(event.target.value as GroupByField)}
          >
            {GROUP_BY_OPTIONS.map((option) => (
              <MenuItem key={option.value} value={option.value}>
                {option.label}
              </MenuItem>
            ))}
          </Select>
        </FormControl>

        <Button
          size="small"
          variant="outlined"
          startIcon={<ViewColumnIcon />}
          aria-haspopup="menu"
          aria-expanded={Boolean(columnsMenuAnchor)}
          aria-controls={columnsMenuAnchor ? "table-columns-menu" : undefined}
          onClick={(event: MouseEvent<HTMLElement>) => setColumnsMenuAnchor(event.currentTarget)}
        >
          Columns
        </Button>
        <Menu
          id="table-columns-menu"
          anchorEl={columnsMenuAnchor}
          open={Boolean(columnsMenuAnchor)}
          onClose={() => setColumnsMenuAnchor(null)}
          anchorOrigin={{ vertical: "bottom", horizontal: "left" }}
          transformOrigin={{ vertical: "top", horizontal: "left" }}
        >
          {TABLE_COLUMNS.map((col) => (
            <MenuItem key={col.id} dense disabled={col.required} onClick={() => toggleColumn(col.id)}>
              <ListItemIcon sx={{ minWidth: 36 }}>
                <Checkbox
                  edge="start"
                  size="small"
                  tabIndex={-1}
                  disableRipple
                  checked={columns.includes(col.id)}
                  disabled={col.required}
                />
              </ListItemIcon>
              <ListItemText primary={col.label} />
            </MenuItem>
          ))}
        </Menu>

        {topGroups.length ? (
          <Stack direction="row" flexWrap="wrap" useFlexGap gap={1}>
            {topGroups.map((group) => (
              <EllipsisChip
                key={group.key}
                clickable
                color={activeChip === group.key ? "primary" : "default"}
                label={`${group.label} (${group.hits.length})`}
                maxWidth={200}
                onClick={() => focusGroup(group.key)}
              />
            ))}
          </Stack>
        ) : null}
      </Stack>

      {groups.map((group) => (
        <Stack key={group.key} id={groupDomId(group.key)} spacing={1}>
          {groupBy !== "none" ? (
            <Typography variant="subtitle1" fontWeight={700}>
              {group.label}{" "}
              <Typography component="span" color="text.secondary">
                ({group.hits.length})
              </Typography>
            </Typography>
          ) : null}
          <TableContainer sx={{ border: 1, borderColor: "divider", borderRadius: 2, overflowX: "auto", bgcolor: TABLE_BG }}>
            <Table size="small" sx={{ tableLayout: "fixed", width: tableWidth, minWidth: "100%" }}>
              <TableHead>
                <TableRow>
                  {visibleColumns.map((col) => (
                    <TableCell
                      key={col.id}
                      align={col.id === "cells" ? "right" : "left"}
                      sx={{
                        width: col.width,
                        verticalAlign: "bottom",
                        fontWeight: 600,
                        px: 1,
                        ...(col.id === "title" ? { ...STICKY_SX, zIndex: 2 } : {}),
                      }}
                    >
                      {col.label}
                    </TableCell>
                  ))}
                </TableRow>
              </TableHead>
              <TableBody>
                {group.hits.map((hit) => (
                  <HitRow key={`${group.key}-${hit.objectID}`} hit={hit} columns={visibleColumns} />
                ))}
              </TableBody>
            </Table>
          </TableContainer>
        </Stack>
      ))}
    </Stack>
  );
}
