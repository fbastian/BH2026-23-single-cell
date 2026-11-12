import { FIELD_TITLES, termLabel, type FacetField } from "../catalogQuery";
import type { ExperimentHit } from "./Hit";

export type GroupByField = "none" | FacetField;

export const GROUP_BY_OPTIONS: { value: GroupByField; label: string }[] = [
  { value: "none", label: "None" },
  ...(
    [
      "organism",
      "tissue",
      "cell_type",
      "development_stage",
      "assay",
      "disease",
      "suspension_type",
      "source",
      "publication",
    ] as FacetField[]
  ).map((field) => ({ value: field as GroupByField, label: FIELD_TITLES[field] })),
];

export type HitGroup = {
  key: string;
  label: string;
  hits: ExperimentHit[];
};

/** Group current page hits under section keys. Multi-valued fields explode into each value. */
export function groupHits(hits: ExperimentHit[], field: GroupByField): HitGroup[] {
  if (field === "none") {
    return [{ key: "__all__", label: "All results", hits }];
  }

  const buckets = new Map<string, ExperimentHit[]>();
  for (const hit of hits) {
    for (const value of new Set(hit.values[field])) {
      const list = buckets.get(value);
      if (list) list.push(hit);
      else buckets.set(value, [hit]);
    }
  }

  return [...buckets.entries()]
    .map(([key, grouped]) => ({ key, label: termLabel(field, key), hits: grouped }))
    .sort((a, b) => a.label.localeCompare(b.label, undefined, { sensitivity: "base" }));
}
