import { ToggleButton, ToggleButtonGroup } from "@mui/material";
import type { CountUnit } from "../catalogQuery";

const COUNT_STORAGE_KEY = "wp1-explore-pivot-count";

/** Initial count unit: URL (?count=publications), else the last choice, else datasets. */
export function loadCountBy(): CountUnit {
  const fromUrl = new URLSearchParams(window.location.search).get("count");
  if (fromUrl === "publications" || fromUrl === "datasets") return fromUrl;
  try {
    return localStorage.getItem(COUNT_STORAGE_KEY) === "publications" ? "publications" : "datasets";
  } catch {
    return "datasets";
  }
}

export function saveCountBy(countBy: CountUnit): void {
  try {
    localStorage.setItem(COUNT_STORAGE_KEY, countBy);
  } catch {
    // ignore
  }
}

/**
 * What the numbers count, in the side filters, the pivot and the result summary: datasets, or
 * publications (the datasets of a publication count once).
 */
export function CountToggle({ value, onChange }: { value: CountUnit; onChange: (value: CountUnit) => void }) {
  return (
    <ToggleButtonGroup
      exclusive
      size="small"
      color="primary"
      value={value}
      aria-label="Count"
      onChange={(_, next: CountUnit | null) => {
        if (next) onChange(next);
      }}
    >
      <ToggleButton value="datasets" title="Count datasets (side filters, pivot)">
        Datasets
      </ToggleButton>
      <ToggleButton
        value="publications"
        title="Count publications (side filters, pivot): the datasets of a publication (same DOI, or same CELLxGENE collection) count once; a dataset without publication information counts as one study"
      >
        Publications
      </ToggleButton>
    </ToggleButtonGroup>
  );
}
