import { useEffect, useMemo, useRef, useState } from "react";
import {
  Alert,
  AppBar,
  Button,
  CircularProgress,
  Container,
  Divider,
  Drawer,
  IconButton,
  List,
  ListItem,
  Pagination,
  Snackbar,
  Stack,
  TextField,
  Toolbar,
  Tooltip,
  Typography,
} from "@mui/material";
import { styled } from "@mui/material/styles";
import ChevronLeftIcon from "@mui/icons-material/ChevronLeft";
import FilterListIcon from "@mui/icons-material/FilterList";
import SearchIcon from "@mui/icons-material/Search";
import LinkIcon from "@mui/icons-material/Link";
import InputAdornment from "@mui/material/InputAdornment";
import { EllipsisChip } from "./components/EllipsisChip";
import { FacetList } from "./components/FacetList";
import { ResultsPanel } from "./components/ResultsPanel";
import { ViewToggle, type ResultsView } from "./components/ViewToggle";
import { CountToggle, loadCountBy, saveCountBy } from "./components/CountToggle";
import { ExportButton } from "./components/ExportButton";
import { setUrlParams, urlList, urlParam } from "./urlState";
import {
  FACET_FIELDS,
  FIELD_TITLES,
  initCatalog,
  searchRows,
  termLabel,
  type ExperimentRow,
  type FacetField,
  type CountUnit,
  type FacetSelection,
} from "./catalogQuery";

const TITLE = "Single-cell datasets across species";

const PAGE_SIZE = 24;
const DRAWER_WIDTH = 300;

const Main = styled("main", {
  shouldForwardProp: (prop) => prop !== "open",
})<{ open?: boolean }>(({ theme, open }) => ({
  flexGrow: 1,
  minWidth: 0,
  transition: theme.transitions.create("margin", {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
  marginLeft: `-${DRAWER_WIDTH}px`,
  ...(open && {
    transition: theme.transitions.create("margin", {
      easing: theme.transitions.easing.easeOut,
      duration: theme.transitions.duration.enteringScreen,
    }),
    marginLeft: 0,
  }),
}));

const DrawerHeader = styled(Toolbar)(({ theme }) => ({
  minHeight: 48,
  paddingLeft: theme.spacing(1),
  paddingRight: theme.spacing(1),
  justifyContent: "flex-end",
  gap: theme.spacing(0.5),
}));

/** Filter panels open at start; the others are collapsed. */
const OPEN_PANELS = new Set<FacetField>([
  "organism",
  "tissue",
  "cell_type",
  "development_stage",
  "suspension_type",
]);

function toggleValue(values: string[] | undefined, value: string): string[] {
  const current = values ?? [];
  return current.includes(value)
    ? current.filter((item) => item !== value)
    : [...current, value];
}

function facetsFromUrl(): FacetSelection {
  const out: FacetSelection = {};
  for (const field of FACET_FIELDS) {
    const values = urlList(field);
    if (values.length) out[field] = values;
  }
  return out;
}

/** Pivot view by default. */
function viewFromUrl(): ResultsView {
  const v = urlParam("view");
  return v === "cards" || v === "table" ? v : "pivot";
}

function pageFromUrl(): number {
  return Math.max(0, (Number(urlParam("page")) || 1) - 1);
}

export default function App() {
  const [ready, setReady] = useState(false);
  const [bootError, setBootError] = useState<string | null>(null);
  // The initial state comes from the URL (see urlState.ts), which is kept up to date below, so
  // that the page can be shared as a link: ?view=table&q=...&tissue=UBERON:0000955,...&page=2
  const [query, setQuery] = useState(() => urlParam("q") ?? "");
  const [debouncedQuery, setDebouncedQuery] = useState(() => urlParam("q") ?? "");
  const [facets, setFacets] = useState<FacetSelection>(facetsFromUrl);
  const [view, setView] = useState<ResultsView>(viewFromUrl);
  const [page, setPage] = useState(pageFromUrl);
  const [linkCopied, setLinkCopied] = useState(false);
  // Incremented by the back/forward buttons (remounts the results, which re-read the URL) and
  // by "Clear" (resets the search boxes and trees of the side filters).
  const [historyVersion, setHistoryVersion] = useState(0);
  const [clearVersion, setClearVersion] = useState(0);
  const [rows, setRows] = useState<ExperimentRow[]>([]);
  const [found, setFound] = useState(0);
  const [foundPublications, setFoundPublications] = useState(0);
  // What the numbers count (side filters, pivot, summary): datasets or publications.
  const [countBy, setCountBy] = useState<CountUnit>(loadCountBy);
  useEffect(() => saveCountBy(countBy), [countBy]);
  const [searchError, setSearchError] = useState<string | null>(null);
  const [searching, setSearching] = useState(false);
  const [filtersOpen, setFiltersOpen] = useState(true);

  useEffect(() => {
    let cancelled = false;
    initCatalog()
      .then(() => {
        if (!cancelled) {
          setReady(true);
          setBootError(null);
        }
      })
      .catch((err: unknown) => {
        if (!cancelled) {
          setBootError(err instanceof Error ? err.message : String(err));
        }
      });
    return () => {
      cancelled = true;
    };
  }, []);

  useEffect(() => {
    const timer = window.setTimeout(() => setDebouncedQuery(query), 180);
    return () => window.clearTimeout(timer);
  }, [query]);

  // Back to the first page when the search or the filters change (not at start: the page can
  // come from the URL). Compared by value, robust to effects run twice in development.
  const searchKey = JSON.stringify([debouncedQuery, facets]);
  const previousSearchKey = useRef(searchKey);
  useEffect(() => {
    if (previousSearchKey.current === searchKey) return;
    previousSearchKey.current = searchKey;
    setPage(0);
  }, [searchKey]);

  // The state in the URL, for sharing links and the back button: a history entry per change,
  // except at start (normalising the URL) and for the search text alone (typing).
  const lastWritten = useRef<string | null>(null);
  useEffect(() => {
    const written = JSON.stringify([view, countBy, page, facets]);
    const mode = lastWritten.current === null || lastWritten.current === written ? "replace" : "push";
    lastWritten.current = written;
    setUrlParams(
      {
        view,
        count: countBy,
        q: debouncedQuery.trim() || null,
        page: page > 0 ? String(page + 1) : null,
        ...Object.fromEntries(FACET_FIELDS.map((field) => [field, facets[field] ?? null])),
      },
      mode
    );
  }, [view, countBy, debouncedQuery, page, facets]);

  // Back/forward buttons: the state comes back from the URL.
  useEffect(() => {
    const onPopState = () => {
      const q = urlParam("q") ?? "";
      const nextFacets = facetsFromUrl();
      const nextView = viewFromUrl();
      const nextPage = pageFromUrl();
      const nextCount = loadCountBy();
      previousSearchKey.current = JSON.stringify([q, nextFacets]);
      lastWritten.current = JSON.stringify([nextView, nextCount, nextPage, nextFacets]);
      setQuery(q);
      setDebouncedQuery(q);
      setFacets(nextFacets);
      setView(nextView);
      setPage(nextPage);
      setCountBy(nextCount);
      setHistoryVersion((v) => v + 1);
    };
    window.addEventListener("popstate", onPopState);
    return () => window.removeEventListener("popstate", onPopState);
  }, []);

  const copyLink = () => {
    navigator.clipboard?.writeText(window.location.href).then(() => setLinkCopied(true), () => undefined);
  };

  useEffect(() => {
    if (!ready || view === "pivot") return;
    let cancelled = false;
    setSearching(true);
    searchRows(debouncedQuery, facets, {
      offset: page * PAGE_SIZE,
      limit: PAGE_SIZE,
    })
      .then((result) => {
        if (!cancelled) {
          setRows(result.rows);
          setFound(result.found);
          setFoundPublications(result.publications);
          setSearchError(null);
          setSearching(false);
        }
      })
      .catch((err: unknown) => {
        if (!cancelled) {
          setSearchError(err instanceof Error ? err.message : String(err));
          setSearching(false);
        }
      });
    return () => {
      cancelled = true;
    };
  }, [ready, debouncedQuery, facets, page, view]);

  const activeRefinements = useMemo(
    () =>
      (Object.entries(facets) as [FacetField, string[]][]).flatMap(([field, values]) =>
        values.map((value) => ({ field, value }))
      ),
    [facets]
  );

  const pageCount = Math.max(1, Math.ceil(found / PAGE_SIZE));

  const onToggleFacet = (field: FacetField, value: string) => {
    setFacets((prev) => {
      const next = { ...prev, [field]: toggleValue(prev[field], value) };
      if (!next[field]?.length) {
        delete next[field];
      }
      return next;
    });
  };

  /**
   * Restrict a field to values (pivot cell or column, table group; several for a merged column
   * of homologous terms). It replaces the field's selection instead of adding to it: values of a
   * field combine with OR, so adding 'forebrain' to 'brain' would still list all brain datasets.
   * The values are within the current selection (the pivot axes and the groups only show
   * selected values and their descendants).
   */
  const onApplyFacet = (field: FacetField, values: string[]) => {
    if (values.length) setFacets((prev) => ({ ...prev, [field]: values }));
  };

  if (bootError) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Typography variant="h1" fontSize={{ xs: "2rem", md: "2.75rem" }} gutterBottom>
          {TITLE}
        </Typography>
        <Alert severity="error">
          {bootError}
          <Typography component="span" display="block" mt={1}>
            Export the Parquet files first: <code>python scripts/export_explore.py</code> (see WP1/README.md)
          </Typography>
        </Alert>
      </Container>
    );
  }

  if (!ready) {
    return (
      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Typography variant="h1" fontSize={{ xs: "2rem", md: "2.75rem" }} gutterBottom>
          {TITLE}
        </Typography>
        <Stack direction="row" spacing={1.5} alignItems="center">
          <CircularProgress size={22} />
          <Typography color="text.secondary">Loading catalog…</Typography>
        </Stack>
      </Container>
    );
  }

  return (
    <Stack sx={{ minHeight: "100vh" }}>
      <AppBar
        position="relative"
        color="transparent"
        elevation={0}
        sx={{
          zIndex: (theme) => theme.zIndex.drawer + 1,
          borderBottom: 1,
          borderColor: "divider",
          bgcolor: "rgba(255,255,255,0.72)",
          backdropFilter: "blur(10px)",
          color: "text.primary",
        }}
      >
        <Container maxWidth="lg" sx={{ py: 2.5 }}>
          <Typography variant="h1" fontSize={{ xs: "1.85rem", md: "2.5rem" }} gutterBottom>
            {TITLE}
          </Typography>
          <Typography color="text.secondary" maxWidth="42rem">
            Number of datasets per species and tissue, cell type or developmental stage, from
            CELLxGENE, scFAIR resources and Fenix, with ontology roll-up (composite-metazoan
            Uberon). Filters on tissue, cell type and stage include more specific terms.
          </Typography>
        </Container>
      </AppBar>

      <Stack direction="row" sx={{ flex: 1, overflow: "hidden" }}>
        <Drawer
          variant="persistent"
          anchor="left"
          open={filtersOpen}
          sx={{
            width: DRAWER_WIDTH,
            flexShrink: 0,
            "& .MuiDrawer-paper": {
              position: "relative",
              width: DRAWER_WIDTH,
              boxSizing: "border-box",
              border: 0,
              borderRight: 1,
              borderColor: "divider",
              bgcolor: "rgba(255,255,255,0.88)",
            },
          }}
        >
          <DrawerHeader variant="dense">
            <Typography variant="subtitle2" fontWeight={700} sx={{ flex: 1, pl: 0.5 }}>
              Filters
              {activeRefinements.length ? ` (${activeRefinements.length})` : ""}
            </Typography>
            <Button
              size="small"
              disabled={!activeRefinements.length && !query}
              onClick={() => {
                setQuery("");
                setFacets({});
                setClearVersion((v) => v + 1);
              }}
            >
              Clear
            </Button>
            <Tooltip title="Collapse filters">
              <IconButton
                size="small"
                edge="end"
                aria-label="Collapse filters"
                onClick={() => setFiltersOpen(false)}
              >
                <ChevronLeftIcon />
              </IconButton>
            </Tooltip>
          </DrawerHeader>
          <Divider />
          <List dense disablePadding sx={{ overflowY: "auto", px: 1.5, py: 1 }}>
            {FACET_FIELDS.map((field) => (
              <ListItem key={field} disableGutters disablePadding sx={{ display: "block" }}>
                <FacetList
                  field={field}
                  title={FIELD_TITLES[field]}
                  query={debouncedQuery}
                  facets={facets}
                  selected={facets[field] ?? []}
                  open={OPEN_PANELS.has(field)}
                  onToggle={(value) => onToggleFacet(field, value)}
                  countBy={countBy}
                  resetSignal={clearVersion}
                />
              </ListItem>
            ))}
          </List>
        </Drawer>

        <Main open={filtersOpen}>
          <Container maxWidth="xl" sx={{ py: 2.5, pb: 8 }}>
            <Stack
              direction={{ xs: "column", sm: "row" }}
              spacing={1.5}
              alignItems={{ sm: "center" }}
              mb={2}
            >
              {!filtersOpen ? (
                <Button
                  variant="outlined"
                  size="small"
                  startIcon={<FilterListIcon />}
                  onClick={() => setFiltersOpen(true)}
                  aria-expanded={false}
                  color={activeRefinements.length ? "primary" : "inherit"}
                  sx={{ flexShrink: 0 }}
                >
                  Filters
                  {activeRefinements.length ? ` (${activeRefinements.length})` : ""}
                </Button>
              ) : null}
              <TextField
                fullWidth
                size="small"
                type="search"
                placeholder="Search titles, publications, DOIs, species, tissues, cell types…"
                value={query}
                onChange={(event) => setQuery(event.target.value)}
                InputProps={{
                  startAdornment: (
                    <InputAdornment position="start">
                      <SearchIcon fontSize="small" />
                    </InputAdornment>
                  ),
                }}
              />
              <Tooltip title="Copy a link to this view (filters, search, pivot, ...): the address of the page is kept up to date">
                <IconButton size="small" aria-label="Copy a link to this view" onClick={copyLink}>
                  <LinkIcon fontSize="small" />
                </IconButton>
              </Tooltip>
              <Snackbar
                open={linkCopied}
                autoHideDuration={2500}
                onClose={() => setLinkCopied(false)}
                message="Link to this view copied to the clipboard"
                anchorOrigin={{ vertical: "bottom", horizontal: "center" }}
              />
              <CountToggle value={countBy} onChange={setCountBy} />
              <ViewToggle value={view} onChange={setView} />
            </Stack>

            {activeRefinements.length ? (
              <Stack direction="row" flexWrap="wrap" useFlexGap gap={1} mb={2}>
                {activeRefinements.map(({ field, value }) => (
                  <EllipsisChip
                    key={`${field}:${value}`}
                    label={`${FIELD_TITLES[field]}: ${termLabel(field, value)}`}
                    maxWidth={220}
                    onDelete={() => onToggleFacet(field, value)}
                  />
                ))}
              </Stack>
            ) : null}

            <Stack direction="row" spacing={2} alignItems="center" mb={1.5}>
              <Typography variant="body2" color="text.secondary">
                {view === "pivot"
                  ? "Pivot view"
                  : searching
                    ? "Searching…"
                    : countBy === "publications"
                      ? `${foundPublications.toLocaleString()} publications (${found.toLocaleString()} datasets)`
                      : `${found.toLocaleString()} datasets`}
              </Typography>
              {view !== "pivot" && found > 0 ? <ExportButton query={debouncedQuery} facets={facets} /> : null}
            </Stack>

            {searchError ? (
              <Alert severity="error" sx={{ mb: 2 }}>
                {searchError}
              </Alert>
            ) : null}

            <ResultsPanel
              key={historyVersion}
              view={view}
              rows={rows}
              query={debouncedQuery}
              facets={facets}
              onOpenHits={() => setView("table")}
              onApplyFacet={onApplyFacet}
              countBy={countBy}
            />

            {view !== "pivot" && found > PAGE_SIZE ? (
              <Stack alignItems="center" mt={3}>
                <Pagination
                  color="primary"
                  page={page + 1}
                  count={pageCount}
                  onChange={(_, next) => setPage(next - 1)}
                />
              </Stack>
            ) : null}
          </Container>
        </Main>
      </Stack>
    </Stack>
  );
}
