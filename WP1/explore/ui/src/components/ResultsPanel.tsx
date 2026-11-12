import { useEffect, useMemo, useState, type MouseEvent, type ReactNode } from "react";
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
  AttributeChips,
  Hit,
  SourceChips,
  formatCompactCount,
  hitTitle,
  publicationLabel,
  publicationTooltip,
  type ExperimentHit,
} from "./Hit";
import { PivotHeatmap } from "./PivotHeatmap";
import type { ResultsView } from "./ViewToggle";

const GROUP_STORAGE_KEY = "wp1-explore-group-by";
// v2: the Publication column was added to the default columns (2026-11-12)
const COLUMNS_STORAGE_KEY = "wp1-explore-table-columns-v2";

type TableColumnId = "title" | "source" | "cells" | FacetField;

const TABLE_CHIP_MAX = 72;

const FIELD_COLUMNS: FacetField[] = [
  "organism",
  "tissue",
  "cell_type",
  "development_stage",
  "assay",
  "disease",
  "suspension_type",
];

const TABLE_COLUMNS: {
  id: TableColumnId;
  label: string;
  /** Fixed width for compact columns; omit so Title can expand. */
  width?: number;
  required?: boolean;
}[] = [
  { id: "title", label: "Title", required: true },
  { id: "source", label: "Resource", width: 96 },
  { id: "publication", label: "Publication", width: 180 },
  ...FIELD_COLUMNS.map((field) => ({ id: field as TableColumnId, label: FIELD_TITLES[field], width: 108 })),
  { id: "cells", label: "Cells", width: 76 },
];

const DEFAULT_COLUMNS: TableColumnId[] = ["title", "source", "publication", "organism", "tissue", "cell_type", "cells"];

function columnCellSx(col: (typeof TABLE_COLUMNS)[number]) {
  return {
    width: col.width ?? "auto",
    maxWidth: col.width ?? "none",
    whiteSpace: "nowrap" as const,
    overflow: "hidden",
    textOverflow: "ellipsis",
    verticalAlign: "middle",
  };
}

function TruncatedCell({
  text,
  width,
  children,
  align,
}: {
  text: string;
  width?: number;
  children?: ReactNode;
  align?: "left" | "right";
}) {
  return (
    <TableCell align={align} sx={{ width: width ?? "auto", maxWidth: width, whiteSpace: "nowrap", verticalAlign: "middle" }}>
      <Tooltip title={text} enterDelay={400}>
        <Typography component="span" display="block" noWrap variant="inherit" sx={{ overflow: "hidden", textOverflow: "ellipsis" }}>
          {children ?? text}
        </Typography>
      </Tooltip>
    </TableCell>
  );
}

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

function ChipCell({ width, children }: { width: number; children: ReactNode }) {
  return (
    <TableCell sx={{ width, maxWidth: width, whiteSpace: "nowrap", verticalAlign: "middle", overflow: "hidden", px: 0.75 }}>
      {children}
    </TableCell>
  );
}

function HitRow({ hit, columns }: { hit: ExperimentHit; columns: typeof TABLE_COLUMNS }) {
  const title = hitTitle(hit);
  return (
    <TableRow hover>
      {columns.map((col) => {
        if (col.id === "title") {
          return (
            <TruncatedCell key={col.id} text={publicationLabel(hit) ? `${title} · ${publicationLabel(hit)}` : title}>
              {hit.explorer_url ? (
                <Link href={hit.explorer_url} target="_blank" rel="noreferrer" underline="hover">
                  {title}
                </Link>
              ) : (
                title
              )}
            </TruncatedCell>
          );
        }
        if (col.id === "source") {
          return (
            <ChipCell key={col.id} width={col.width!}>
              <Stack direction="row" gap={0.5} sx={{ overflow: "hidden" }}>
                <SourceChips hit={hit} />
              </Stack>
            </ChipCell>
          );
        }
        if (col.id === "publication") {
          const label = publicationLabel(hit);
          return (
            <TruncatedCell key={col.id} text={label ? publicationTooltip(hit) : "No publication information"} width={col.width}>
              {label && hit.publication_url ? (
                <Link href={hit.publication_url} target="_blank" rel="noreferrer" underline="hover">
                  {label}
                </Link>
              ) : (
                label ?? "—"
              )}
            </TruncatedCell>
          );
        }
        if (col.id === "cells") {
          return (
            <TruncatedCell key={col.id} text={hit.cell_count.toLocaleString()} width={col.width} align="right">
              {formatCompactCount(hit.cell_count)}
            </TruncatedCell>
          );
        }
        const field = col.id as FacetField;
        return (
          <ChipCell key={col.id} width={col.width!}>
            <AttributeChips values={labelsOf(field, hit.values[field])} chipMaxWidth={TABLE_CHIP_MAX} wrap={false} />
          </ChipCell>
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
          <TableContainer sx={{ border: 1, borderColor: "divider", borderRadius: 2 }}>
            <Table size="small" sx={{ tableLayout: "fixed", width: "100%" }}>
              <TableHead>
                <TableRow>
                  {visibleColumns.map((col) => (
                    <TableCell key={col.id} align={col.id === "cells" ? "right" : "left"} sx={columnCellSx(col)}>
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
