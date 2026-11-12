import { useEffect, useMemo, useState } from "react";
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Breadcrumbs,
  Button,
  Checkbox,
  IconButton,
  Link,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import {
  allValues,
  childrenOf,
  facetCounts,
  isHierarchical,
  topTermsOf,
  visibleValues,
  termLabel,
  type CountUnit,
  type FacetField,
  type FacetSelection,
  type HierarchicalField,
} from "../catalogQuery";

type Props = {
  field: FacetField;
  title: string;
  query: string;
  facets: FacetSelection;
  selected: string[];
  open?: boolean;
  limit?: number;
  onToggle: (value: string) => void;
  countBy: CountUnit;
};

/** Values to list: by count, then label. */
function byCount(counts: Map<string, number>, field: FacetField) {
  return (a: string, b: string) =>
    (counts.get(b) ?? 0) - (counts.get(a) ?? 0) ||
    termLabel(field, a).localeCompare(termLabel(field, b));
}

/**
 * Facet filter. Flat fields list their values. Hierarchical fields (tissue, cell type,
 * stage) show the children of the current term, with counts including descendants; the arrow
 * drills down, the breadcrumb goes back up, and the search box searches all terms. Ticking a
 * term matches the datasets annotated with it or with any of its descendants.
 */
export function FacetList({
  field,
  title,
  query,
  facets,
  selected,
  open = true,
  limit = 12,
  onToggle,
  countBy,
}: Props) {
  const hierarchical = isHierarchical(field);
  const [counts, setCounts] = useState<Map<string, number>>(new Map());
  const [filter, setFilter] = useState("");
  const [expanded, setExpanded] = useState(open);
  const [showMore, setShowMore] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // Drill-down path in the hierarchy, starting at the top term if there is one visible root
  // (empty = list the top terms: the roots, or the visible children of hidden roots).
  const [topPath] = useState<string[]>(() => {
    const tops = hierarchical ? topTermsOf(field as HierarchicalField) : [];
    return tops.length === 1 ? tops : [];
  });
  const [path, setPath] = useState<string[]>(topPath);

  useEffect(() => {
    if (!expanded) return;
    let cancelled = false;
    const handle = window.setTimeout(() => {
      facetCounts(field, query, facets, countBy)
        .then((result) => {
          if (!cancelled) {
            setCounts(result);
            setError(null);
          }
        })
        .catch((err: unknown) => {
          if (!cancelled) setError(err instanceof Error ? err.message : String(err));
        });
    }, 120);
    return () => {
      cancelled = true;
      window.clearTimeout(handle);
    };
  }, [field, query, facets, expanded, countBy]);

  const needle = filter.trim().toLowerCase();
  const values = useMemo(() => {
    const present = (v: string) => (counts.get(v) ?? 0) > 0;
    let list: string[];
    if (needle) {
      list = visibleValues(field).filter(
        (v) => present(v) && termLabel(field, v).toLowerCase().includes(needle)
      );
    } else if (hierarchical) {
      const node = path[path.length - 1];
      list = (node ? childrenOf(field as HierarchicalField, node) : topTermsOf(field as HierarchicalField))
        .filter(present);
    } else {
      list = [...counts.keys()].filter(present);
    }
    return list.sort(byCount(counts, field));
  }, [counts, needle, hierarchical, path, field]);

  const selectedSet = new Set(selected);
  const visible = showMore ? values : values.slice(0, limit);
  const node = path[path.length - 1];
  const searchable = hierarchical || allValues(field).length > limit;

  const hasChildren = (v: string) =>
    hierarchical &&
    childrenOf(field as HierarchicalField, v).some((c) => (counts.get(c) ?? 0) > 0);

  const row = (value: string, bold = false) => {
    const labelId = `facet-${field}-${value}`;
    return (
      <ListItem
        key={value}
        disablePadding
        dense
        secondaryAction={
          !needle && hasChildren(value) && value !== node ? (
            <Tooltip title="Show more specific terms">
              <IconButton
                edge="end"
                size="small"
                aria-label={`Show terms under ${termLabel(field, value)}`}
                onClick={() => {
                  setPath((p) => [...p, value]);
                  setShowMore(false);
                }}
              >
                <ChevronRightIcon fontSize="small" />
              </IconButton>
            </Tooltip>
          ) : null
        }
      >
        <ListItemButton
          role={undefined}
          dense
          onClick={() => onToggle(value)}
          sx={{ py: 0.25, pr: hierarchical ? 5 : 1, alignItems: "center" }}
        >
          <ListItemIcon sx={{ minWidth: 36 }}>
            <Checkbox
              edge="start"
              size="small"
              tabIndex={-1}
              disableRipple
              checked={selectedSet.has(value)}
              inputProps={{ "aria-labelledby": labelId }}
              sx={{ p: 0.5 }}
            />
          </ListItemIcon>
          <ListItemText
            id={labelId}
            primary={termLabel(field, value)}
            primaryTypographyProps={{
              variant: "body2",
              fontWeight: bold ? 700 : undefined,
              sx: { wordBreak: "break-word" },
            }}
          />
          <Typography variant="caption" color="text.secondary" sx={{ pl: 1, flexShrink: 0 }}>
            {(counts.get(value) ?? 0).toLocaleString()}
          </Typography>
        </ListItemButton>
      </ListItem>
    );
  };

  return (
    <Accordion
      disableGutters
      elevation={0}
      expanded={expanded}
      onChange={(_, next) => setExpanded(next)}
      sx={{
        bgcolor: "transparent",
        borderBottom: 1,
        borderColor: "divider",
        "&:before": { display: "none" },
      }}
    >
      <AccordionSummary expandIcon={<ExpandMoreIcon />}>
        <Typography variant="subtitle2" fontWeight={700}>
          {title}
          {selected.length ? ` (${selected.length})` : ""}
        </Typography>
      </AccordionSummary>
      <AccordionDetails sx={{ pt: 0, px: 0 }}>
        {searchable ? (
          <TextField
            size="small"
            fullWidth
            type="search"
            placeholder={`Search ${title.toLowerCase()}…`}
            value={filter}
            onChange={(event) => setFilter(event.target.value)}
            sx={{ mb: 1, px: 1.5 }}
          />
        ) : null}
        {hierarchical && !needle && path.length > topPath.length ? (
          <Breadcrumbs maxItems={4} sx={{ px: 1.5, mb: 0.5, fontSize: "0.75rem" }}>
            {topPath.length ? null : (
              <Link component="button" variant="caption" onClick={() => setPath([])}>
                top
              </Link>
            )}
            {path.slice(0, -1).map((id, i) => (
              <Link
                key={id}
                component="button"
                variant="caption"
                onClick={() => setPath((p) => p.slice(0, i + 1))}
              >
                {termLabel(field, id)}
              </Link>
            ))}
            <Typography variant="caption" color="text.primary">
              {termLabel(field, node!)}
            </Typography>
          </Breadcrumbs>
        ) : null}
        {error ? (
          <Typography color="error" variant="caption" display="block" sx={{ mb: 1, px: 1.5 }}>
            {error}
          </Typography>
        ) : null}
        <List dense disablePadding>
          {hierarchical && !needle && node ? row(node, true) : null}
          <List dense disablePadding sx={hierarchical && !needle && node ? { pl: 2 } : undefined}>
            {visible.map((value) => row(value))}
          </List>
        </List>
        {!values.length ? (
          <Typography variant="caption" color="text.secondary" sx={{ px: 1.5 }}>
            No matching values
          </Typography>
        ) : null}
        {values.length > limit ? (
          <Button size="small" onClick={() => setShowMore((v) => !v)} sx={{ mt: 1, ml: 1.5 }}>
            {showMore ? "Show less" : `Show all ${values.length}`}
          </Button>
        ) : null}
      </AccordionDetails>
    </Accordion>
  );
}
