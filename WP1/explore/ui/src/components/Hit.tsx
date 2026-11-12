import { Fragment, useState, type MouseEvent } from "react";
import { Box, Card, CardContent, Chip, Link, Stack, Tooltip, Typography } from "@mui/material";
import ArticleOutlinedIcon from "@mui/icons-material/ArticleOutlined";
import {
  FIELD_TITLES,
  isPublication,
  labelsOf,
  publicationSize,
  termLabel,
  type ExperimentRow,
  type FacetField,
} from "../catalogQuery";

export type ExperimentHit = ExperimentRow;

export function hitTitle(hit: ExperimentHit): string {
  return hit.title || hit.id;
}

/**
 * Compact counts for chips/table.
 * Scaled coefficient under 10 → one decimal (1.5k, 1.3M); 10+ → integer (45k, 34M).
 */
export function formatCompactCount(value: number): string {
  const n = Number.isFinite(value) ? value : 0;
  const abs = Math.abs(n);
  const formatScaled = (scaled: number) => {
    const rounded = Math.round(scaled * 10) / 10;
    return Math.abs(rounded) < 10 ? rounded.toFixed(1) : String(Math.round(scaled));
  };
  if (abs < 1000) return String(Math.round(n));
  if (abs < 1_000_000) return `${formatScaled(n / 1000)}k`;
  if (abs < 1_000_000_000) return `${formatScaled(n / 1_000_000)}M`;
  return `${formatScaled(n / 1_000_000_000)}B`;
}

/**
 * Chip style whose label wraps onto several lines instead of being cut: a chip's label is on
 * one line by default and gets an ellipsis when the chip is wider than its column.
 */
export const WRAP_CHIP_SX = {
  height: "auto",
  maxWidth: "100%",
  "& .MuiChip-label": {
    display: "block",
    whiteSpace: "normal",
    overflowWrap: "break-word",
    py: "2px",
    lineHeight: 1.35,
  },
} as const;

/**
 * A list of values as chips ("bubbles") whose labels wrap: the first `limit` values, then a
 * "+N more" chip that shows the whole list ("show less" folds it again).
 */
export function ValueChips({ values, limit = 6 }: { values: string[]; limit?: number }) {
  const [expanded, setExpanded] = useState(false);
  if (!values.length) {
    return (
      <Typography variant="body2" color="text.secondary">
        —
      </Typography>
    );
  }
  const shown = expanded ? values : values.slice(0, limit);
  const hidden = values.length - shown.length;
  const toggle = (open: boolean) => (event: MouseEvent) => {
    event.stopPropagation();
    setExpanded(open);
  };
  return (
    <Stack direction="row" flexWrap="wrap" useFlexGap gap={0.5}>
      {shown.map((value, i) => (
        <Chip key={`${i}-${value}`} size="small" variant="outlined" label={value} sx={WRAP_CHIP_SX} />
      ))}
      {hidden > 0 ? (
        <Chip size="small" color="primary" variant="outlined" clickable label={`+${hidden} more`} onClick={toggle(true)} />
      ) : null}
      {expanded && values.length > limit ? (
        <Chip size="small" clickable label="show less" onClick={toggle(false)} />
      ) : null}
    </Stack>
  );
}

/** Flags of a dataset (merged across resources, several species) as warning chips. */
export function FlagChips({ hit }: { hit: ExperimentHit }) {
  return (
    <>
      {hit.values.flags.map((flag) => (
        <Chip key={flag} size="small" color="warning" variant="outlined" label={termLabel("flags", flag)} sx={WRAP_CHIP_SX} />
      ))}
    </>
  );
}

/** Chips of the resources of a dataset, each linking to the dataset in that resource. */
export function SourceChips({ hit }: { hit: ExperimentHit }) {
  return (
    <>
      {hit.values.source.map((slug, i) => {
        const url = hit.source_explorer_urls[i];
        return (
          <Chip
            key={slug}
            size="small"
            color="primary"
            label={termLabel("source", slug)}
            sx={WRAP_CHIP_SX}
            {...(url ? { component: "a", href: url, target: "_blank", rel: "noreferrer", clickable: true } : {})}
          />
        );
      })}
    </>
  );
}

/** Label of the dataset's publication, or null without publication information. */
export function publicationLabel(hit: ExperimentHit): string | null {
  const pub = hit.values.publication[0];
  return isPublication(pub) ? termLabel("publication", pub) : null;
}

/** Tooltip of a publication: its label and its number of datasets in the explorer. */
export function publicationTooltip(hit: ExperimentHit): string {
  const pub = hit.values.publication[0];
  const n = isPublication(pub) ? publicationSize(pub) : 0;
  return `Publication: ${publicationLabel(hit)}${n > 1 ? ` · ${n} datasets in this explorer` : ""}`;
}

/** Chip of the dataset's publication (link to the DOI or collection); nothing without one. */
export function PublicationChip({ hit }: { hit: ExperimentHit }) {
  const label = publicationLabel(hit);
  if (!label) return null;
  const url = hit.publication_url;
  return (
    <Tooltip title={publicationTooltip(hit)} enterDelay={300}>
      <Chip
        size="small"
        variant="outlined"
        icon={<ArticleOutlinedIcon />}
        label={label}
        sx={WRAP_CHIP_SX}
        {...(url ? { component: "a", href: url, target: "_blank", rel: "noreferrer", clickable: true } : {})}
      />
    </Tooltip>
  );
}

/** Fields shown on a card, as on the side filters (resource, publication and flags are chips). */
export const DETAIL_FIELDS: FacetField[] = [
  "organism",
  "tissue",
  "cell_type",
  "development_stage",
  "assay",
  "disease",
  "sex",
  "suspension_type",
  "tissue_type",
];

export function Hit({ hit }: { hit: ExperimentHit }) {
  const title = hitTitle(hit);
  return (
    <Card variant="outlined" sx={{ height: "100%" }}>
      <CardContent sx={{ p: 1.5, "&:last-child": { pb: 1.5 } }}>
        <Typography variant="subtitle1" component="h3" fontWeight={700} lineHeight={1.25} sx={{ mb: 0.75 }}>
          {hit.explorer_url ? (
            <Link href={hit.explorer_url} target="_blank" rel="noreferrer" underline="hover">
              {title}
            </Link>
          ) : (
            title
          )}
        </Typography>
        <Stack direction="row" flexWrap="wrap" useFlexGap gap={0.5} mb={1.25}>
          <SourceChips hit={hit} />
          <Tooltip title={hit.cell_count.toLocaleString()} enterDelay={400}>
            <Chip size="small" label={`${formatCompactCount(hit.cell_count)} cells`} />
          </Tooltip>
          <FlagChips hit={hit} />
          <PublicationChip hit={hit} />
          {hit.h5ad_url ? (
            <Chip size="small" variant="outlined" component="a" href={hit.h5ad_url} clickable label="H5AD" />
          ) : null}
        </Stack>
        <Box
          sx={{
            display: "grid",
            gridTemplateColumns: "max-content minmax(0, 1fr)",
            columnGap: 1.5,
            rowGap: 0.75,
            alignItems: "start",
          }}
        >
          {DETAIL_FIELDS.map((field) => (
            <Fragment key={field}>
              <Typography variant="caption" color="text.secondary" fontWeight={600} sx={{ pt: "3px" }}>
                {FIELD_TITLES[field]}
              </Typography>
              <ValueChips values={labelsOf(field, hit.values[field])} />
            </Fragment>
          ))}
        </Box>
      </CardContent>
    </Card>
  );
}
