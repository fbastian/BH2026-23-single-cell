import type { ReactNode } from "react";
import { Card, CardContent, Chip, Link, Stack, Tooltip, Typography } from "@mui/material";
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
import { EllipsisChip } from "./EllipsisChip";

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

/** Attribute lists: at most `limit` value chips, then a single +n chip listing the rest. */
export function AttributeChips({
  values,
  limit = 3,
  chipMaxWidth = 110,
  wrap = true,
}: {
  values?: string[];
  limit?: number;
  chipMaxWidth?: number;
  wrap?: boolean;
}) {
  if (!values?.length) {
    return (
      <Typography variant="caption" color="text.secondary">
        —
      </Typography>
    );
  }
  const shown = values.slice(0, limit);
  const rest = values.slice(limit);
  return (
    <Stack
      direction="row"
      flexWrap={wrap ? "wrap" : "nowrap"}
      useFlexGap
      gap={0.5}
      sx={wrap ? undefined : { overflow: "hidden" }}
    >
      {shown.map((value) => (
        <EllipsisChip key={value} size="small" label={value} variant="outlined" maxWidth={chipMaxWidth} />
      ))}
      {rest.length ? (
        <Tooltip title={rest.join(" · ")}>
          <Chip size="small" label={`+${rest.length}`} variant="outlined" />
        </Tooltip>
      ) : null}
    </Stack>
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
export function PublicationChip({ hit, maxWidth = 260 }: { hit: ExperimentHit; maxWidth?: number }) {
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
        sx={{ maxWidth }}
        {...(url ? { component: "a", href: url, target: "_blank", rel: "noreferrer", clickable: true } : {})}
      />
    </Tooltip>
  );
}

function FacetBlock({ label, children }: { label: string; children: ReactNode }) {
  return (
    <Stack spacing={0.25} minWidth={0}>
      <Typography variant="caption" color="text.secondary" fontWeight={600} lineHeight={1.2}>
        {label}
      </Typography>
      {children}
    </Stack>
  );
}

const CARD_FIELDS: FacetField[] = [
  "organism",
  "tissue",
  "cell_type",
  "development_stage",
  "assay",
  "disease",
  "suspension_type",
  "sex",
];

export function Hit({ hit }: { hit: ExperimentHit }) {
  const title = hitTitle(hit);
  return (
    <Card variant="outlined" sx={{ height: "100%" }}>
      <CardContent sx={{ p: 1.25, "&:last-child": { pb: 1.25 } }}>
        <Typography
          variant="subtitle1"
          component="h3"
          fontWeight={700}
          lineHeight={1.25}
          sx={{
            display: "-webkit-box",
            WebkitLineClamp: 2,
            WebkitBoxOrient: "vertical",
            overflow: "hidden",
            mb: 0.75,
          }}
        >
          {hit.explorer_url ? (
            <Link href={hit.explorer_url} target="_blank" rel="noreferrer" underline="hover">
              {title}
            </Link>
          ) : (
            title
          )}
        </Typography>
        <Stack direction="row" flexWrap="wrap" useFlexGap gap={0.5} mb={1}>
          <SourceChips hit={hit} />
          <Tooltip title={hit.cell_count.toLocaleString()} enterDelay={400}>
            <Chip size="small" label={`${formatCompactCount(hit.cell_count)} cells`} />
          </Tooltip>
          {hit.values.flags.map((flag) => (
            <Chip key={flag} size="small" color="warning" variant="outlined" label={termLabel("flags", flag)} />
          ))}
          <PublicationChip hit={hit} />
          {hit.h5ad_url ? (
            <Chip size="small" variant="outlined" component="a" href={hit.h5ad_url} clickable label="H5AD" />
          ) : null}
        </Stack>
        <Stack
          direction="row"
          flexWrap="wrap"
          useFlexGap
          gap={1}
          sx={{ "& > *": { flex: "1 1 calc(50% - 8px)", minWidth: 0, maxWidth: "100%" } }}
        >
          {CARD_FIELDS.map((field) => (
            <FacetBlock key={field} label={FIELD_TITLES[field]}>
              <AttributeChips values={labelsOf(field, hit.values[field])} />
            </FacetBlock>
          ))}
        </Stack>
      </CardContent>
    </Card>
  );
}
