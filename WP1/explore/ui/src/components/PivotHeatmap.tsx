import { useEffect, useMemo, useRef, useState } from "react";
import { setUrlParams, urlList, urlParam } from "../urlState";
import { HeatMap } from "@nivo/heatmap";
import type { AxisTickProps } from "@nivo/axes";
import {
  Alert,
  Autocomplete,
  Breadcrumbs,
  CircularProgress,
  FormControl,
  InputLabel,
  Link,
  MenuItem,
  Paper,
  Select,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";
import {
  visibleValues,
  pathTo,
  stripExact,
  SELECTED_NODE,
  EXACT_PREFIX,
  HOMOLOGY_PREFIX,
  filterValuesFor,
  hasValue,
  homologyGroup,
  withHomologs,
  FACET_FIELDS,
  FIELD_TITLES,
  SELF_PREFIX,
  childrenOf,
  isHierarchical,
  pivotCounts,
  type CountUnit,
  rootsOf,
  termLabel,
  type FacetField,
  type FacetSelection,
  type HierarchicalField,
  type PivotAxis,
  type PivotResult,
} from "../catalogQuery";

const STORAGE_KEY = "wp1-explore-pivot-axes";
const TOP_OPTIONS = [15, 30, 60, 120];
const MAX_LABEL = 34;

/** Axis state: field; for hierarchical fields, roll-up on/off and the drill-down path. */
type AxisState = { field: FacetField; rollup: boolean; path: string[] };

const NIVO_THEME = {
  text: { fontFamily: '"IBM Plex Sans", "Segoe UI", sans-serif', fontSize: 11, fill: "#14212b" },
  axis: {
    ticks: {
      text: { fontFamily: '"IBM Plex Sans", "Segoe UI", sans-serif', fontSize: 11, fill: "#5b6b78" },
    },
  },
  tooltip: { container: { fontFamily: '"IBM Plex Sans", "Segoe UI", sans-serif', fontSize: 12 } },
} as const;

const LABEL_ON_LIGHT = "#14212b";
const LABEL_ON_DARK = "#f4fffb";

/** Pick light/dark label text from cell fill for readable contrast. */
function contrastingLabelColor(background: string): string {
  let r = 0;
  let g = 0;
  let b = 0;
  const rgb = background.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/i);
  if (rgb) {
    [r, g, b] = [Number(rgb[1]), Number(rgb[2]), Number(rgb[3])];
  } else {
    const hex = background.replace("#", "").trim();
    const full = hex.length === 3 ? hex.split("").map((ch) => ch + ch).join("") : hex;
    if (full.length >= 6) {
      r = Number.parseInt(full.slice(0, 2), 16);
      g = Number.parseInt(full.slice(2, 4), 16);
      b = Number.parseInt(full.slice(4, 6), 16);
    }
  }
  const channel = (c: number) => {
    const s = c / 255;
    return s <= 0.03928 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  };
  const luminance = 0.2126 * channel(r) + 0.7152 * channel(g) + 0.0722 * channel(b);
  return luminance > 0.42 ? LABEL_ON_LIGHT : LABEL_ON_DARK;
}

/**
 * Where the roll-up of a field starts, given the filter selection on that field: the root
 * without selection, the selected term if one, the selected terms (SELECTED_NODE) if several.
 */
function startPath(field: FacetField, scope: string[]): string[] {
  if (!isHierarchical(field)) return [];
  const ids = [...new Set(scope.map(stripExact))];
  if (!ids.length) return rootsOf(field).slice(0, 1);
  return ids.length === 1 ? pathTo(field, ids[0]) : [SELECTED_NODE];
}

/** path [] = start from the filter selection (see startPath). */
function newAxis(field: FacetField, rollup: boolean): AxisState {
  return { field, rollup: rollup && isHierarchical(field), path: [] };
}

/** An axis from the URL: rows=organism, cols=tissue&cols_mode=rollup&cols_path=UBERON:...,... */
function axisFromUrl(prefix: "rows" | "cols"): AxisState | null {
  const field = urlParam(prefix);
  if (!field || !(FACET_FIELDS as readonly string[]).includes(field)) return null;
  const f = field as FacetField;
  const path = urlList(`${prefix}_path`).filter((id) => id === SELECTED_NODE || hasValue(f, id));
  return { field: f, rollup: urlParam(`${prefix}_mode`) !== "raw" && isHierarchical(f), path };
}

/** Initial axes: from the URL, else the last choice (localStorage), else species x tissue roll-up. */
function loadAxes(): { row: AxisState; col: AxisState } {
  const stored = loadStoredAxes();
  return { row: axisFromUrl("rows") ?? stored.row, col: axisFromUrl("cols") ?? stored.col };
}

function loadStoredAxes(): { row: AxisState; col: AxisState } {
  try {
    const parsed = JSON.parse(localStorage.getItem(STORAGE_KEY) ?? "null") as {
      row?: { field?: string; rollup?: boolean };
      col?: { field?: string; rollup?: boolean };
    } | null;
    const ok = (f?: string): f is FacetField => (FACET_FIELDS as readonly string[]).includes(f ?? "");
    if (parsed && ok(parsed.row?.field) && ok(parsed.col?.field)) {
      return {
        row: newAxis(parsed.row!.field as FacetField, Boolean(parsed.row!.rollup)),
        col: newAxis(parsed.col!.field as FacetField, Boolean(parsed.col!.rollup)),
      };
    }
  } catch {
    // ignore
  }
  return { row: newAxis("organism", false), col: newAxis("tissue", true) };
}

function shorten(text: string): string {
  return text.length > MAX_LABEL ? `${text.slice(0, MAX_LABEL - 1)}…` : text;
}

/** Tooltip note of the column named after the roll-up node itself. */
const SELF_NOTE = " annotated with this term itself and with no term below it";
const SELF_NOTE_MERGED = " annotated with one of these terms itself and with no term below it";

const CELL_WIDTH = 46;
const CELL_HEIGHT = 28;
const CHART_SIDE = 12;

/** Row labels as HTML, pinned on the left while the heatmap scrolls horizontally. */
function RowLabels({
  rows,
  totals,
  field,
  width,
  top,
  height,
  canDrill,
  onDrill,
  unit,
}: {
  rows: string[];
  totals: number[];
  field: FacetField;
  width: number;
  top: number;
  height: number;
  canDrill: (value: string) => boolean;
  onDrill: (value: string) => void;
  unit: string;
}) {
  return (
    <Stack
      sx={{
        position: "sticky",
        left: 0,
        zIndex: 1,
        flexShrink: 0,
        width,
        height,
        pt: `${top}px`,
        bgcolor: "rgba(255,255,255,0.96)",
        borderRight: 1,
        borderColor: "divider",
      }}
    >
      {rows.map((value, i) => {
        const label = termLabel(field, value);
        const clickable = canDrill(value);
        const id = value.startsWith(SELF_PREFIX) ? value.slice(SELF_PREFIX.length) : value;
        return (
          <Typography
            key={value}
            variant="caption"
            noWrap
            title={`${label} (${id}): ${(totals[i] ?? 0).toLocaleString()} ${unit}${value.startsWith(SELF_PREFIX) ? SELF_NOTE : ""}${clickable ? " · click to show only the terms under it" : ""}`}
            onClick={clickable ? () => onDrill(value) : undefined}
            sx={{
              height: CELL_HEIGHT,
              lineHeight: `${CELL_HEIGHT}px`,
              textAlign: "right",
              pr: 1,
              pl: 1,
              color: clickable ? "primary.main" : "text.secondary",
              textDecoration: clickable ? "underline" : undefined,
              cursor: clickable ? "pointer" : "default",
            }}
          >
            {shorten(label)}
          </Typography>
        );
      })}
    </Stack>
  );
}

type Props = {
  query: string;
  facets: FacetSelection;
  /** What the cells count (set by the Count switch of the page). */
  countBy: CountUnit;
  /** Restrict fields to values (pivot cell or column): the datasets listed are those counted. */
  onApplyCell: (cells: { field: FacetField; values: string[] }[]) => void;
};

function AxisControls({
  id,
  title,
  axis,
  onChange,
}: {
  id: string;
  title: string;
  axis: AxisState;
  onChange: (axis: AxisState) => void;
}) {
  return (
    <Stack direction="row" spacing={1} alignItems="center">
      <FormControl size="small" sx={{ minWidth: 170 }}>
        <InputLabel id={`${id}-label`}>{title}</InputLabel>
        <Select
          labelId={`${id}-label`}
          label={title}
          value={axis.field}
          onChange={(event) => onChange(newAxis(event.target.value as FacetField, axis.rollup))}
        >
          {FACET_FIELDS.map((field) => (
            <MenuItem key={field} value={field}>
              {FIELD_TITLES[field]}
            </MenuItem>
          ))}
        </Select>
      </FormControl>
      {isHierarchical(axis.field) ? (
        <ToggleButtonGroup
          exclusive
          size="small"
          color="primary"
          value={axis.rollup ? "rollup" : "raw"}
          onChange={(_, next: string | null) => {
            if (next) onChange({ ...axis, rollup: next === "rollup", path: [] });
          }}
        >
          <ToggleButton value="raw" title="Terms as annotated (after remapping)">
            Raw
          </ToggleButton>
          <ToggleButton value="rollup" title="Ontology roll-up: all terms under the current term, each counting the datasets annotated with it or with a more specific term">
            Roll-up
          </ToggleButton>
        </ToggleButtonGroup>
      ) : null}
    </Stack>
  );
}

function AxisBreadcrumbs({
  title,
  axis,
  path,
  scope,
  onChange,
}: {
  title: string;
  axis: AxisState;
  /** Current roll-up path (the axis path, or the start from the filter selection). */
  path: string[];
  /** Filter selection on the field (named in the breadcrumb when several terms are selected). */
  scope: string[];
  onChange: (axis: AxisState) => void;
}) {
  const field = axis.field;
  const crumb = (id: string) =>
    id === SELECTED_NODE ? scope.map((s) => termLabel(field, s)).join(" + ") : termLabel(field, id);
  const options = useMemo(
    () =>
      isHierarchical(field)
        ? visibleValues(field)
            .filter((id) => childrenOf(field, id).length > 0)
            .sort((a, b) => termLabel(field, a).localeCompare(termLabel(field, b)))
        : [],
    [field]
  );
  if (!axis.rollup || !isHierarchical(field)) return null;
  return (
    <Stack direction="row" spacing={1} alignItems="center" flexWrap="wrap" useFlexGap>
      <Typography variant="body2" color="text.secondary" sx={{ flexShrink: 0 }}>
        {title}: under
      </Typography>
      <Breadcrumbs maxItems={6} sx={{ fontSize: "0.875rem" }}>
        {path.map((id, i) =>
          i === path.length - 1 ? (
            <Typography key={id} variant="body2" color="text.primary" fontWeight={600}>
              {crumb(id)}
            </Typography>
          ) : (
            <Link
              key={id}
              component="button"
              variant="body2"
              onClick={() => onChange({ ...axis, path: path.slice(0, i + 1) })}
            >
              {crumb(id)}
            </Link>
          )
        )}
      </Breadcrumbs>
      <Autocomplete
        size="small"
        sx={{ width: 280 }}
        options={options}
        value={null}
        blurOnSelect
        clearOnBlur
        getOptionLabel={(id) => termLabel(field, id)}
        renderOption={(props, id) => (
          <li {...props} key={id}>
            <Typography variant="body2">
              {termLabel(field, id)}{" "}
              <Typography component="span" variant="caption" color="text.secondary">
                {id}
              </Typography>
            </Typography>
          </li>
        )}
        onChange={(_, id) => {
          if (id) onChange({ ...axis, path: pathTo(field, id) });
        }}
        renderInput={(params) => <TextField {...params} placeholder="Jump to a term…" />}
      />
    </Stack>
  );
}

export function PivotHeatmap({ query, facets, countBy, onApplyCell }: Props) {
  const [axes, setAxes] = useState(loadAxes);
  const [top, setTop] = useState(() => {
    const n = Number(urlParam("top"));
    return TOP_OPTIONS.includes(n) ? n : 30;
  });
  const unit = countBy === "publications" ? "publications" : "datasets";

  const [data, setData] = useState<PivotResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // The selection of the column field also brings the terms homologous to the selected terms
  // (pivot only: side filters and lists keep the literal selection), so that their merged
  // column appears when both have data with the other filters.
  const homologs = useMemo(
    () => withHomologs(axes.col.field, [...new Set(facets[axes.col.field] ?? [])]),
    [facets, axes.col.field]
  );
  const pivotFacets = useMemo(
    () => (homologs.added.length ? { ...facets, [axes.col.field]: homologs.values } : facets),
    [facets, homologs, axes.col.field]
  );

  // Filter selection on the field of each axis: the axis shows only these values (and their
  // descendants), and the roll-up starts from them. A new selection resets the drill-down.
  const rowScope = useMemo(() => [...new Set(pivotFacets[axes.row.field] ?? [])], [pivotFacets, axes.row.field]);
  const colScope = useMemo(() => [...new Set(pivotFacets[axes.col.field] ?? [])], [pivotFacets, axes.col.field]);
  const rowScopeKey = rowScope.join("\n");
  const colScopeKey = colScope.join("\n");
  // Only when the selection changes, not at start (the drill-down can come from the URL).
  const previousRowScope = useRef(rowScopeKey);
  const previousColScope = useRef(colScopeKey);
  useEffect(() => {
    if (previousRowScope.current === rowScopeKey) return;
    previousRowScope.current = rowScopeKey;
    setAxes((prev) => (prev.row.path.length ? { ...prev, row: { ...prev.row, path: [] } } : prev));
  }, [rowScopeKey]);
  useEffect(() => {
    if (previousColScope.current === colScopeKey) return;
    previousColScope.current = colScopeKey;
    setAxes((prev) => (prev.col.path.length ? { ...prev, col: { ...prev.col, path: [] } } : prev));
  }, [colScopeKey]);

  const rowPath = axes.row.path.length ? axes.row.path : startPath(axes.row.field, rowScope);
  const colPath = axes.col.path.length ? axes.col.path : startPath(axes.col.field, colScope);
  const pathOf = (which: "row" | "col") => (which === "row" ? rowPath : colPath);
  const rowNode = axes.row.rollup ? rowPath[rowPath.length - 1] ?? null : null;
  const colNode = axes.col.rollup ? colPath[colPath.length - 1] ?? null : null;
  const rowAxis: PivotAxis = { field: axes.row.field, node: rowNode, scope: rowScope };
  const colAxis: PivotAxis = { field: axes.col.field, node: colNode, scope: colScope };

  useEffect(() => {
    try {
      const keep = (a: AxisState) => ({ field: a.field, rollup: a.rollup });
      localStorage.setItem(STORAGE_KEY, JSON.stringify({ row: keep(axes.row), col: keep(axes.col) }));
    } catch {
      // ignore
    }
  }, [axes]);

  // The pivot state in the URL, for sharing links and the back button (a history entry per
  // change; the first write of the component only normalises the URL).
  const pivotWritten = useRef(false);
  useEffect(() => {
    const mode = (a: AxisState) => (isHierarchical(a.field) ? (a.rollup ? "rollup" : "raw") : null);
    const historyMode = pivotWritten.current ? "push" : "replace";
    pivotWritten.current = true;
    setUrlParams({
      rows: axes.row.field,
      rows_mode: mode(axes.row),
      rows_path: axes.row.path,
      cols: axes.col.field,
      cols_mode: mode(axes.col),
      cols_path: axes.col.path,
      top: String(top),
    }, historyMode);
  }, [axes, top]);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);
    const timer = window.setTimeout(() => {
      pivotCounts(
        rowAxis,
        colAxis,
        query,
        pivotFacets,
        { maxRows: top, maxCols: top, countBy }
      )
        .then((result) => {
          if (!cancelled) {
            setData(result);
            setLoading(false);
          }
        })
        .catch((err: unknown) => {
          if (!cancelled) {
            setError(err instanceof Error ? err.message : String(err));
            setLoading(false);
          }
        });
    }, 150);
    return () => {
      cancelled = true;
      window.clearTimeout(timer);
    };
  }, [axes.row.field, axes.col.field, rowNode, colNode, rowScope, colScope, query, pivotFacets, top, countBy]);

  const chartData = useMemo(() => {
    if (!data) return [];
    return data.rows.map((row, rowIndex) => ({
      id: row,
      data: data.cols.map((col, colIndex) => ({ x: col, y: data.matrix[rowIndex]?.[colIndex] ?? 0 })),
    }));
  }, [data]);

  /** A value can be drilled into if its axis is rolled up and it has children. */
  const canDrill = (axis: AxisState, value: string) =>
    axis.rollup &&
    isHierarchical(axis.field) &&
    !value.startsWith(SELF_PREFIX) &&
    !value.startsWith(EXACT_PREFIX) &&
    childrenOf(axis.field as HierarchicalField, value).length > 0;

  const drill = (which: "row" | "col", value: string) =>
    setAxes((prev) => ({ ...prev, [which]: { ...prev[which], path: [...pathOf(which), value] } }));

  const tickRenderer = (which: "row" | "col") => {
    const axis = axes[which];
    const totals = data ? (which === "row" ? data.rowTotals : data.colTotals) : [];
    const order = data ? (which === "row" ? data.rows : data.cols) : [];
    return (tick: AxisTickProps<string>) => {
      const value = String(tick.value);
      const label = termLabel(axis.field, value);
      const clickable = canDrill(axis, value);
      const total = totals[order.indexOf(value)] ?? 0;
      const base = value.startsWith(SELF_PREFIX) ? value.slice(SELF_PREFIX.length) : value;
      const group = base.startsWith(HOMOLOGY_PREFIX) ? homologyGroup(axis.field, base) : undefined;
      const selfNote = !value.startsWith(SELF_PREFIX) ? "" : group ? SELF_NOTE_MERGED : SELF_NOTE;
      const title = group
        ? `${label}: homologous terms merged into one column (homologous in ${group.taxa}; Bgee anatomical similarity annotations): ${total.toLocaleString()} ${unit}${selfNote}`
        : `${label} (${base}): ${total.toLocaleString()} ${unit}${selfNote}`;
      return (
        <g transform={`translate(${tick.x},${tick.y})`}>
          <line x1={0} x2={tick.lineX} y1={0} y2={tick.lineY} stroke="#b8c4cc" />
          <text
            transform={`translate(${tick.textX},${tick.textY}) rotate(${tick.rotate ?? 0})`}
            textAnchor={tick.textAnchor as "start" | "middle" | "end"}
            dominantBaseline={tick.textBaseline as "central"}
            onClick={clickable ? () => drill(which, value) : undefined}
            style={{
              fontFamily: '"IBM Plex Sans", "Segoe UI", sans-serif',
              fontSize: 11,
              fill: clickable ? "#0f6a5c" : "#5b6b78",
              textDecoration: clickable ? "underline" : undefined,
              cursor: clickable ? "pointer" : "default",
              fontStyle: group ? "italic" : undefined,
            }}
          >
            <title>
              {title}
              {clickable ? " · click to show only the terms under it" : ""}
            </title>
            {shorten(label)}
          </text>
        </g>
      );
    };
  };

  /** Button under a column: lists all the datasets of the column (table view). */
  const columnButton = (tick: AxisTickProps<string>) => {
    const value = String(tick.value);
    const label = termLabel(axes.col.field, value);
    const total = data ? data.colTotals[data.cols.indexOf(value)] ?? 0 : 0;
    return (
      <g
        transform={`translate(${tick.x},${tick.y + 8})`}
        onClick={() => onApplyCell([{ field: axes.col.field, values: filterValuesFor(colAxis, value) }])}
        style={{ cursor: "pointer" }}
        role="button"
        aria-label={`List the datasets of ${label}`}
      >
        <title>{`List all the datasets of ${label} (${total.toLocaleString()} ${unit})`}</title>
        <rect x={-18} y={0} width={36} height={18} rx={9} fill="#ffffff" stroke="#0f6a5c" />
        <text
          x={0}
          y={9}
          textAnchor="middle"
          dominantBaseline="central"
          style={{ fontFamily: '"IBM Plex Sans", "Segoe UI", sans-serif', fontSize: 10, fill: "#0f6a5c" }}
        >
          list
        </text>
      </g>
    );
  };

  const longest = (values: string[], field: FacetField) =>
    Math.max(4, ...values.map((v) => Math.min(MAX_LABEL, termLabel(field, v).length)));
  // Fixed cell size: the table keeps its size and scrolls horizontally; row labels are an
  // HTML column pinned on the left (sticky), aligned with the rows of the heatmap.
  const labelWidth = data ? Math.min(260, 24 + longest(data.rows, axes.row.field) * 6.6) : 140;
  const marginTop = data ? Math.min(200, 30 + longest(data.cols, axes.col.field) * 4.6) : 110;
  // Below the cells: one "list" button per column, then the colour legend.
  const marginBottom = 86;
  const chartWidth = data ? Math.max(300, data.cols.length * CELL_WIDTH + 2 * CHART_SIDE) : 300;
  const chartHeight = data ? data.rows.length * CELL_HEIGHT + marginTop + marginBottom : 360;

  const setAxis = (which: "row" | "col") => (axis: AxisState) =>
    setAxes((prev) => ({ ...prev, [which]: axis }));

  return (
    <Stack spacing={1.5}>
      <Stack direction={{ xs: "column", md: "row" }} spacing={2} alignItems={{ md: "center" }} flexWrap="wrap" useFlexGap>
        <AxisControls id="pivot-rows" title="Rows" axis={axes.row} onChange={setAxis("row")} />
        <AxisControls id="pivot-cols" title="Columns" axis={axes.col} onChange={setAxis("col")} />
        <FormControl size="small" sx={{ minWidth: 110 }}>
          <InputLabel id="pivot-top-label">Show top</InputLabel>
          <Select
            labelId="pivot-top-label"
            label="Show top"
            value={top}
            onChange={(event) => setTop(Number(event.target.value))}
          >
            {TOP_OPTIONS.map((n) => (
              <MenuItem key={n} value={n}>
                {n}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
      </Stack>
      <AxisBreadcrumbs title="Rows" axis={axes.row} path={rowPath} scope={rowScope} onChange={setAxis("row")} />
      <AxisBreadcrumbs title="Columns" axis={axes.col} path={colPath} scope={colScope} onChange={setAxis("col")} />

      {data ? (
        <Typography variant="body2" color="text.secondary">
          {countBy === "publications"
            ? `${data.units.toLocaleString()} publications (${data.found.toLocaleString()} datasets) · a publication counts once per row, column and cell, whatever its number of datasets; a dataset without publication information counts as one study`
            : `${data.found.toLocaleString()} datasets`}{" "}
          · a dataset counts in every row and column it is annotated with (several species,
          tissues, ...) · click a cell to list its datasets, or "list" under a column to list all the
          datasets of the column
          {axes.row.rollup || axes.col.rollup ? " · roll-up: a dataset also counts for all the more general terms of its annotations; click an underlined label to show only the terms under it" : ""}
        </Typography>
      ) : null}
      {homologs.added.length ? (
        <Alert severity="info" sx={{ py: 0 }}>
          The columns also include the terms homologous to your selection (Bgee anatomical
          similarity annotations):{" "}
          {homologs.added
            .map(([m, s]) => `${termLabel(axes.col.field, m)} (homologous to ${termLabel(axes.col.field, s)})`)
            .join(", ")}
          . Their merged column appears when they have data with the other filters.
        </Alert>
      ) : null}

      {loading ? (
        <Stack direction="row" spacing={1} alignItems="center">
          <CircularProgress size={18} />
          <Typography variant="body2" color="text.secondary">
            Computing pivot…
          </Typography>
        </Stack>
      ) : null}
      {error ? <Alert severity="error">{error}</Alert> : null}

      {data && !loading && !data.rows.length ? (
        <Alert severity="info">No dataset matches.</Alert>
      ) : null}

      {data && !loading && data.rows.length ? (
        <Paper variant="outlined" sx={{ width: "100%", overflowX: "auto", bgcolor: "rgba(255,255,255,0.72)" }}>
          <Stack direction="row" sx={{ width: labelWidth + chartWidth }}>
            <RowLabels
              rows={data.rows}
              totals={data.rowTotals}
              field={axes.row.field}
              width={labelWidth}
              top={marginTop}
              height={chartHeight}
              canDrill={(value) => canDrill(axes.row, value)}
              onDrill={(value) => drill("row", value)}
              unit={unit}
            />
            <HeatMap
              width={chartWidth}
              height={chartHeight}
              data={chartData}
              margin={{ top: marginTop, right: CHART_SIDE, bottom: marginBottom, left: CHART_SIDE }}
              theme={NIVO_THEME}
              forceSquare={false}
              xInnerPadding={0.08}
              yInnerPadding={0.08}
              axisTop={{ tickSize: 4, tickPadding: 6, tickRotation: -45, renderTick: tickRenderer("col") }}
              axisRight={null}
              axisBottom={{ tickSize: 0, tickPadding: 0, renderTick: columnButton }}
              axisLeft={null}
              colors={{
                type: "sequential",
                colors: ["#e8f2ef", "#0f6a5c"],
                minValue: 0,
                maxValue: Math.max(1, data.maxCell),
              }}
              emptyColor="#f3f6f4"
              borderColor={{ from: "color", modifiers: [["darker", 0.15]] }}
              borderWidth={1}
              borderRadius={2}
              enableLabels
              labelTextColor={(cell) => contrastingLabelColor(cell.color)}
              hoverTarget="cell"
              animate={false}
              isInteractive
              onClick={(cell) => {
                if (!cell.value) return;
                onApplyCell([
                  { field: axes.row.field, values: filterValuesFor(rowAxis, String(cell.serieId)) },
                  { field: axes.col.field, values: filterValuesFor(colAxis, String(cell.data.x)) },
                ]);
              }}
              tooltip={({ cell }) => (
                <Paper elevation={3} sx={{ px: 1.25, py: 0.75 }}>
                  <Typography variant="caption" display="block" color="text.secondary">
                    {termLabel(axes.row.field, String(cell.serieId))} ×{" "}
                    {termLabel(axes.col.field, String(cell.data.x))}
                  </Typography>
                  <Typography variant="body2" fontWeight={700}>
                    {(cell.value ?? 0).toLocaleString()} {unit}
                  </Typography>
                  {countBy === "publications" ? (
                    <Typography variant="caption" display="block" color="text.secondary">
                      {(
                        data.datasetMatrix[data.rows.indexOf(String(cell.serieId))]?.[
                          data.cols.indexOf(String(cell.data.x))
                        ] ?? 0
                      ).toLocaleString()}{" "}
                      datasets
                    </Typography>
                  ) : null}
                </Paper>
              )}
              legends={[
                {
                  anchor: "bottom",
                  translateX: 0,
                  translateY: 58,
                  length: 160,
                  thickness: 8,
                  direction: "row",
                  tickPosition: "after",
                  tickSize: 3,
                  tickSpacing: 4,
                  tickOverlap: false,
                  title: countBy === "publications" ? "Publications" : "Datasets",
                  titleAlign: "start",
                  titleOffset: 4,
                },
              ]}
            />
          </Stack>
        </Paper>
      ) : null}
    </Stack>
  );
}
