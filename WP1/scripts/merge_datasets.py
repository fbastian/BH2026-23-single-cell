"""Merge dataset metadata exported from several sources into one deduplicated JSON.

Inputs are given with --input FORMAT=PATH (repeatable). The order defines the priority: the
first input is preferred for scalar fields and for choosing the primary record of a merge.
Each FORMAT has a reader (see READERS) that turns the parsed JSON file into a list of
(input dataset id, record in the common format below). To support a new source, add a reader
function and register it in READERS.

Common record format (the CELLxGENE/scFAIR format): `dataset_id`, `source` ({slug, name}),
`source_reference_id`, `source_url`, `organism`, `tissue`, `cell_type`, `development_stage`,
`sex`, `disease`, `assay` (lists of {label, ontology_term_id}), `suspension_type` (list of
strings), `cell_count`, `collection_id`, `collection_doi`, `explorer_url`, `assets`, ...

Deduplication
  1. Same record: same (source slug, source_reference_id), e.g. a CELLxGENE dataset present in
     both the CELLxGENE export and the scFAIR export.
  2. Same data in different resources: same DOI, same set of organisms and same cell count,
     different sources. A match is used only if it is unambiguous (each source appears at most
     once in the group); otherwise the records are left separate and counted in the summary.

Merging keeps the most complete information: list fields are the union of the records' values
(deduplicated by ontology_term_id, or label when there is no term id, assets by URL; tissues of
the same term with different tissue_type, e.g. tissue and organoid, are kept apart); other
fields are taken from the highest-priority record that has a value. A record of an outdated
version of the dataset (its `dataset_version_id` differs from that of the highest-priority
record having one) only fills fields that are empty in the current records: put the freshest
source first. scFAIR copies of CELLxGENE datasets get their `dataset_version_id` from the name
of their H5AD file. Missing term fields become empty lists. DOIs are stored without URL
encoding and without doi.org prefix (case kept), and compared case-insensitively. Output
records have
  dataset_id  "<slug>:<source_reference_id>" of the primary (highest-priority) record
  sources     one entry per merged record (source, source_reference_id, source_url,
              explorer_url, input_format, input_dataset_id)
  merged      true when records from several sources were merged (case 2)
Entries of term lists holding several terms joined with ' || ' are split.

Usage: python WP1/scripts/merge_datasets.py \\
           --input cxg=WP1/data/datasets_CxG.json --input scfair=WP1/data/datasets_scFAIR.json \\
           --input fenix=WP1/data/datasets_fenix.json \\
           --out WP1/output/datasets_merged.json [--report WP1/output/merged_datasets.tsv]
"""

import argparse
import collections
import csv
import json
import re
import urllib.parse
from pathlib import Path

TERM_FIELDS = ("organism", "tissue", "cell_type", "development_stage", "sex", "disease", "assay",
               "self_reported_ethnicity")
SOURCE_NAMES = {"cxg": "CELLxGENE", "scp": "Single Cell Portal", "asap": "ASAP", "bgee": "Bgee"}


def split_joined_terms(entries):
    """Split entries like {label: 'a || b', ontology_term_id: 'X:1 || X:2'} into one per term."""
    out = []
    for e in entries:
        ids, labels = str(e["ontology_term_id"]).split(" || "), str(e["label"]).split(" || ")
        if len(ids) > 1 and len(ids) == len(labels):
            out.extend({**e, "label": lab, "ontology_term_id": i} for i, lab in zip(ids, labels))
        else:
            out.append(e)
    return out


def clean_doi(doi):
    """DOI without URL encoding (e.g. '10.1016%2Fj.cell...' from Single Cell Portal), surrounding
    spaces and doi.org/doi: prefix; case kept. None if empty."""
    if not doi:
        return None
    doi = urllib.parse.unquote(str(doi)).strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "http://dx.doi.org/", "doi:"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix):]
    return doi or None


def normalise(rec):
    rec = dict(rec)
    for f in TERM_FIELDS:
        if isinstance(rec.get(f), list):
            rec[f] = split_joined_terms([e for e in rec[f] if isinstance(e, dict)])
        elif rec.get(f) is None:
            rec[f] = []  # missing field (e.g. self_reported_ethnicity outside CELLxGENE)
    if "collection_doi" in rec:
        rec["collection_doi"] = clean_doi(rec["collection_doi"])
    return rec


def read_cxg(data):
    """CELLxGENE curation API export: a list of records."""
    out = []
    for rec in data:
        rec = dict(rec)
        rec["source"] = {"slug": "cxg", "name": SOURCE_NAMES["cxg"]}
        rec["source_reference_id"] = rec["dataset_id"]
        rec["source_url"] = f"https://cellxgene.cziscience.com/collections/{rec['collection_id']}"
        out.append((rec["dataset_id"], rec))
    return out


# CELLxGENE H5AD files are named after the dataset version id
CXG_H5AD_URL = re.compile(r"^https://datasets\.cellxgene\.cziscience\.com/([0-9a-f-]{36})\.h5ad$")


def read_scfair(data):
    """scFAIR portal export: a list of records that already have source and source_reference_id.

    scFAIR copies of CELLxGENE datasets have no dataset_version_id; it is read from the name of
    their H5AD file, so that copies of an older version are recognised as outdated."""
    out = []
    for rec in data:
        rec = dict(rec)
        if rec["source"]["slug"] == "cxg" and not rec.get("dataset_version_id"):
            versions = {m.group(1) for a in rec.get("assets") or []
                        if (m := CXG_H5AD_URL.match(a.get("url") or ""))}
            if len(versions) == 1:
                rec["dataset_version_id"] = versions.pop()
        out.append((rec["dataset_id"], rec))
    return out


# Fields of the common format that the 'fenix' records carry under the same name
FENIX_COPIED_FIELDS = (
    "title", "cell_count", "assets", "explorer_url", "organism", "tissue", "cell_type",
    "development_stage", "sex", "disease", "assay", "suspension_type", "self_reported_ethnicity",
    "donor_id", "is_primary_data", "mean_genes_per_cell", "primary_cell_count", "published_at",
    "revised_at", "schema_version", "spatial", "tombstone", "x_approximate_distribution",
    "batch_condition",
)


def read_fenix(data):
    """CELLxGENE Discover search export ({"data": [...], "total": n}). These are CELLxGENE
    datasets (experiment_id = dataset id, study_id = collection id), so they get source 'cxg'
    and are merged with the other CELLxGENE records of the same dataset."""
    out = []
    for rec in data["data"]:
        new = {f: rec[f] for f in FENIX_COPIED_FIELDS if f in rec}
        new["tissue"] = [{k: v for k, v in t.items() if not k.startswith("top-organ")}
                         for t in new.get("tissue", [])]
        doi = re.search(r"doi\.org/(\S+)", rec.get("citation") or "")
        new["collection_doi"] = doi.group(1).rstrip(".,;)") if doi else None
        new["dataset_id"] = rec["experiment_id"]
        new["dataset_version_id"] = rec.get("experiment_version_id")
        new["collection_id"] = rec["study_id"]
        new["source"] = {"slug": "cxg", "name": SOURCE_NAMES["cxg"]}
        new["source_reference_id"] = rec["experiment_id"]
        new["source_url"] = f"https://cellxgene.cziscience.com/collections/{rec['study_id']}"
        out.append((rec["_id"], new))
    return out


READERS = {"cxg": read_cxg, "scfair": read_scfair, "fenix": read_fenix}


def norm_doi(doi):
    """DOI for comparisons: cleaned (clean_doi) and lower-cased (DOIs are case-insensitive)."""
    doi = clean_doi(doi)
    return doi.lower() if doi else None


def is_empty(v):
    return v is None or v == [] or v == {} or v == ""


def item_key(item):
    if isinstance(item, dict) and "ontology_term_id" in item:
        # some entries have no term id (e.g. Bgee assays): fall back to the label
        return item["ontology_term_id"] or f"label:{item.get('label')}"
    if isinstance(item, dict) and "url" in item:
        return item["url"]
    return json.dumps(item, sort_keys=True)


def union_items(lists):
    """Union of list values, deduplicated by item_key, in priority order. CELLxGENE tissue
    entries of the same term with different tissue_type (e.g. 'tissue' and 'organoid') are kept
    apart; an entry without tissue_type (scFAIR) is a duplicate of any entry of the same term."""
    out, positions = [], collections.defaultdict(list)  # key -> indexes in out
    for v in lists:
        for item in v:
            k = item_key(item)
            tt = item.get("tissue_type") if isinstance(item, dict) else None
            types = {out[i].get("tissue_type") if isinstance(out[i], dict) else None
                     for i in positions[k]}
            if positions[k] and (tt is None or tt in types):
                continue
            if types == {None}:  # replace the untyped entry by the typed one
                out[positions[k][0]] = item
                continue
            positions[k].append(len(out))
            out.append(item)
    return out


def merge_records(records):
    """records: list of (priority, input_format, input_dataset_id, normalised record).

    Records of an outdated version of the dataset (a `dataset_version_id` different from the one
    of the highest-priority record that has one) only fill fields that no current record has;
    they never add to lists, so stale terms and files do not leak into the merge.
    """
    records = sorted(records, key=lambda r: r[0])
    primary = records[0][3]
    ref_version = next((r[3]["dataset_version_id"] for r in records
                        if r[3].get("dataset_version_id")), None)

    def outdated(r):
        v = r[3].get("dataset_version_id")
        return bool(v and ref_version and v != ref_version)

    current = [r for r in records if not outdated(r)]
    out = {}
    fields = list(dict.fromkeys(k for r in records for k in r[3]))
    for f in fields:
        values = [r[3][f] for r in current if f in r[3] and not is_empty(r[3][f])]
        if not values:
            values = [r[3][f] for r in records if f in r[3] and not is_empty(r[3][f])][:1]
        if not values:
            out[f] = records[0][3].get(f)
        elif isinstance(values[0], list):
            out[f] = union_items(values)
        else:
            out[f] = values[0]
    slugs = {r[3]["source"]["slug"] for r in records}
    out["dataset_id"] = f"{primary['source']['slug']}:{primary['source_reference_id']}"
    out["source"] = primary["source"]
    out["source_reference_id"] = primary["source_reference_id"]
    out["sources"] = [
        {
            "source": r[3]["source"],
            "source_reference_id": r[3]["source_reference_id"],
            "source_url": r[3].get("source_url"),
            "explorer_url": r[3].get("explorer_url"),
            "input_format": r[1],
            "input_dataset_id": r[2],
        }
        for r in records
    ]
    out["merged"] = len(slugs) > 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--input", action="append", required=True, metavar="FORMAT=PATH",
                    help=f"input file and its format ({', '.join(READERS)}); repeat; first = highest priority")
    ap.add_argument("--out", type=Path, required=True, help="merged JSON to write")
    ap.add_argument("--report", type=Path, help="optional TSV listing cross-resource merges")
    args = ap.parse_args()

    records = []  # (priority, format, original dataset_id, record)
    for priority, spec in enumerate(args.input):
        fmt, _, path = spec.partition("=")
        if fmt not in READERS or not path:
            ap.error(f"--input must be FORMAT=PATH with FORMAT in {sorted(READERS)}: {spec}")
        read = READERS[fmt](json.load(open(path, encoding="utf-8")))
        for input_id, rec in read:
            records.append((priority, fmt, input_id, normalise(rec)))
        print(f"{path}: {len(read)} records ({fmt})")

    # 1. same source record
    groups = collections.defaultdict(list)
    for r in records:
        groups[(r[3]["source"]["slug"], r[3]["source_reference_id"])].append(r)
    merged_same = [merge_records(g) for g in groups.values()]
    n_exact = sum(len(g) > 1 for g in groups.values())

    # 2. same data in different resources
    def cross_key(m):
        doi = norm_doi(m.get("collection_doi"))
        orgs = tuple(sorted(o["ontology_term_id"] for o in m.get("organism", [])))
        if not doi or not orgs or m.get("cell_count") in (None, 0):
            return None
        return doi, orgs, m["cell_count"]

    by_key = collections.defaultdict(list)
    for i, m in enumerate(merged_same):
        k = cross_key(m)
        if k:
            by_key[k].append(i)
    consumed, final, report = set(), [], []
    n_cross = n_ambiguous = 0
    for k, idx in by_key.items():
        slugs = [merged_same[i]["source"]["slug"] for i in idx]
        if len(set(slugs)) < 2:
            continue
        if len(slugs) != len(set(slugs)):
            n_ambiguous += 1
            continue
        # rebuild from the original records so that priorities are respected
        members = [r for i in idx for r in groups[(merged_same[i]["source"]["slug"],
                                                   merged_same[i]["source_reference_id"])]]
        m = merge_records(members)
        final.append(m)
        consumed.update(idx)
        n_cross += 1
        report.append({"doi": k[0], "organisms": ";".join(k[1]), "cell_count": k[2],
                       "dataset_id": m["dataset_id"],
                       "merged_records": ";".join(dict.fromkeys(
                           f"{s['source']['slug']}:{s['source_reference_id']}" for s in m["sources"]))})
    final.extend(m for i, m in enumerate(merged_same) if i not in consumed)
    final.sort(key=lambda m: m["dataset_id"])

    args.out.write_text(json.dumps(final, ensure_ascii=False, indent=1), encoding="utf-8")
    if args.report:
        with args.report.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["doi", "organisms", "cell_count", "dataset_id",
                                              "merged_records"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(sorted(report, key=lambda r: r["dataset_id"]))

    print(f"{len(records)} input records -> {len(final)} datasets in {args.out}")
    print(f"  same-source duplicates merged: {n_exact}")
    print(f"  cross-resource merges (same DOI + organisms + cell count): {n_cross}")
    print(f"  ambiguous cross-resource groups left unmerged: {n_ambiguous}")
    print("  datasets per source:", dict(collections.Counter(m["source"]["slug"] for m in final)))


if __name__ == "__main__":
    main()
