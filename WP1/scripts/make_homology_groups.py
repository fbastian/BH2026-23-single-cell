"""Group terms that are homologous across species, from Bgee's anatomical similarity annotations.

Reads the 'ancestral taxa homology annotations' of BgeeDB/anatomical-similarity-annotations
(TSV columns: HOM ID, HOM name, entity, entity name, taxon ID, taxon name, CIO ID, CIO name,
supporting text). Only multi-entity annotations are used: their entity is several ids
separated by '|' (e.g. UBERON:0002048|UBERON:0006860, lung and swim bladder, homologous in
Gnathostomata), so they relate different terms, which the explore page merges into one pivot
column. Single-entity annotations (a structure homologous across the species of a taxon) relate
nothing that composite-metazoan Uberon does not already merge.

Each entity id is mapped to composite-metazoan: kept if it is a term of composite, else replaced
by the composite terms it is an xref of (e.g. ZFA or FBbt ids); ids mapping to nothing are
reported. For each factor, the members kept are the terms of the ancestor TSV of that factor
(annotated terms and their ancestors). Annotations sharing a member are merged (connected
components: pectoral fin|arm and pectoral fin|forelimb give one group with the three terms), so
that a term is in one group at most. Groups with fewer than two members are dropped.

Output TSV, one row per member:
  factor        tissue, cell_type or development_stage
  group_id      the member ids of the group, joined with '|'
  group_label   the member labels, joined with ' | '
  term_id, term_label, annotated   as in the ancestor TSV
  taxa          ancestral taxa in which the annotations of the group hold ('; ')
  annotations   entity column of these annotations ('; ')
  confidence    their CIO confidence ('; ')

The input files are not modified.

Usage: python WP1/scripts/make_homology_groups.py --annotations ancestral_taxa_homology_annotations.tsv
           --composite composite-metazoan.obo --ancestors term_ancestors.tsv --out term_homology.tsv
"""

import argparse
import collections
import csv
from pathlib import Path

from make_term_remappings import FACTORS, parse_obo


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--annotations", type=Path, required=True,
                    help="ancestral_taxa_homology_annotations.tsv (BgeeDB/anatomical-similarity-annotations)")
    ap.add_argument("--composite", type=Path, required=True,
                    help="path to composite-metazoan.obo (reference ontology)")
    ap.add_argument("--ancestors", type=Path, required=True,
                    help="term ancestor TSV written by make_term_ancestors.py")
    ap.add_argument("--out", type=Path, required=True, help="homology group TSV to write")
    args = ap.parse_args()

    composite = parse_obo(args.composite)
    in_composite = {i for i, t in composite.items() if not t["obsolete"] and (t["name"] or t["parents"])}
    xref_of = collections.defaultdict(set)
    for i in in_composite:
        for x in composite[i]["xrefs"]:
            xref_of[x].add(i)

    def to_composite(entity_id):
        return {entity_id} if entity_id in in_composite else xref_of.get(entity_id, set())

    terms = {}  # (factor, term id) -> ancestor TSV row
    with args.ancestors.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            terms[(r["factor"], r["term_id"])] = r

    with args.annotations.open(encoding="utf-8", newline="") as f:
        annotations = [r for r in csv.DictReader(f, delimiter="\t") if "|" in (r.get("entity") or "")]
    unmapped = sorted({e for a in annotations for e in a["entity"].split("|") if not to_composite(e)})

    rows, n_groups = [], collections.Counter()
    for factor in FACTORS:
        # union-find over the members present in this factor
        parent = {}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        used = []  # (annotation, its members in this factor)
        for a in annotations:
            members = sorted({c for e in a["entity"].split("|") for c in to_composite(e)
                              if (factor, c) in terms})
            if len(members) < 2:
                continue
            used.append((a, members))
            for m in members:
                parent.setdefault(m, m)
            for m in members[1:]:
                parent[find(m)] = find(members[0])
        groups = collections.defaultdict(set)
        for m in parent:
            groups[find(m)].add(m)
        for members in groups.values():
            members = sorted(members)
            anns = [a for a, ms in used if set(ms) <= set(members)]
            group_id = "|".join(members)
            label = " | ".join(terms[(factor, m)]["term_label"] for m in members)
            unique = lambda values: "; ".join(dict.fromkeys(values))
            for m in members:
                rows.append({
                    "factor": factor, "group_id": group_id, "group_label": label,
                    "term_id": m, "term_label": terms[(factor, m)]["term_label"],
                    "annotated": terms[(factor, m)]["annotated"],
                    "taxa": unique(a["taxon name"] for a in anns),
                    "annotations": unique(a["entity"] for a in anns),
                    "confidence": unique(a["CIO name"] for a in anns),
                })
            n_groups[factor] += 1

    rows.sort(key=lambda r: (r["factor"], r["group_label"], r["term_label"]))
    fields = ["factor", "group_id", "group_label", "term_id", "term_label", "annotated", "taxa",
              "annotations", "confidence"]
    with args.out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    if unmapped:
        print(f"warning: {len(unmapped)} annotation ids not in composite (ignored): " + ", ".join(unmapped))
    print(f"{len(annotations)} multi-entity annotations; {sum(n_groups.values())} groups, {len(rows)} terms -> {args.out}")
    for factor in FACTORS:
        labels = sorted({r["group_label"] for r in rows if r["factor"] == factor})
        if labels:
            print(f"  {factor}: {n_groups[factor]} groups: " + "; ".join(labels))


if __name__ == "__main__":
    main()
