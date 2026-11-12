"""Export the merged datasets, with remapped terms and the term hierarchy, for the explore UI.

Reads the merged datasets JSON (merge_datasets.py), the term remapping TSV
(make_term_remappings.py), the term ancestor TSV (make_term_ancestors.py) and the term subset
TSV (make_term_subsets.py), and writes two Parquet files loaded by the browser
(WP1/explore/ui):

--out-datasets: one row per dataset. Columns holding facet values are lists of value ids:
  organism           NCBITaxon ids; subspecies are folded into their species (a label of 3 or
                     more words becomes its first two words, with the id of that species in
                     the datasets, e.g. Mus musculus musculus -> Mus musculus)
  tissue, cell_type, development_stage
                     terms after remapping (remap -> targets, exclude -> dropped); an empty
                     annotation is 'na', remapped to the root of the factor by the TSV
  assay, disease, sex
                     ontology term ids ('label:<label>' for entries without id); 'unknown',
                     'na' and empty annotations are all 'na'
  suspension_type    cell, nucleus or na ('na' and empty are both 'na')
  tissue_type        tissue types of the tissue entries (CELLxGENE only), else 'na'
  source             resources of the dataset (several for datasets merged across resources)
  publication        the publication (or study) of the dataset, one per dataset, so that the
                     page can count publications instead of datasets. Key, in this order:
                     doi:<DOI> (lower case); else, if the dataset's collection has a DOI on
                     other datasets, that DOI; else collection:<resource>:<collection id>;
                     else dataset:<dataset id> (no publication information: the dataset counts
                     as one study). Collection ids are used only for the resources of
                     COLLECTION_RESOURCES (CELLxGENE: in the scFAIR export, the collection id of
                     Bgee and Single Cell Portal datasets is the dataset's own accession, and
                     ASAP has a placeholder, ASAP000000, for all projects), and a collection
                     whose datasets have several DOIs is ignored. Label: DOI label (e.g. 'Siletti et al. (2023) Science'), else
                     collection name, else DOI; the dataset title for dataset: keys; the DOI or
                     id is appended to labels shared by different publications.
  flags              'merged' (merged across resources), 'multi_species' (several species)
  plus: id, title (falls back to collection name, then DOI label or resource name with the
  resource id), cell_count, explorer_url, h5ad_url, doi, doi_label, collection_name,
  publication_url (doi.org link, or collection page; empty for dataset: keys),
  source_explorer_urls (explorer URL per entry of `source`), search_blob (text search).

--out-terms: one row per value of each field: field, id, label, parent_ids, ancestor_ids,
  subsets. For tissue, cell type and developmental stage, these come from the ancestor TSV
  (annotated terms and their ancestors within the factor domain); the UI computes the roll-up
  from them. subsets: the subsets of the term in the subset TSV (e.g. upper_level); the UI
  hides these terms in the filter trees and the roll-up columns. Other fields have no
  hierarchy and no subsets.

The input files are not modified.

Usage: python WP1/scripts/export_explore.py --datasets datasets_merged.json
           --remappings term_remappings.tsv --ancestors term_ancestors.tsv
           --subsets term_subsets.tsv --out-datasets datasets.parquet --out-terms terms.parquet
"""

import argparse
import collections
import csv
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from make_term_ancestors import read_remappings, remapped_terms
from make_term_remappings import FACTORS
from merge_datasets import clean_doi, norm_doi

NA = "na"
NA_LABEL = "NA (not annotated)"
UNKNOWN_IDS = {"unknown", "na"}
TERM_FIELDS = ("assay", "disease", "sex")
FLAG_LABELS = {"merged": "merged across resources", "multi_species": "several species"}
# Resources whose collection_id groups the datasets of a publication (see 'publication' above)
COLLECTION_RESOURCES = ("cxg",)
LIST_COLUMNS = ("organism", *FACTORS, *TERM_FIELDS, "suspension_type", "tissue_type",
                "source", "publication", "flags", "source_explorer_urls")


def read_ancestors(path):
    """{factor: {term_id: (label, [parent ids], [ancestor ids])}} from make_term_ancestors.py."""
    terms = collections.defaultdict(dict)
    with path.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            split = lambda s: s.split(";") if s else []
            terms[r["factor"]][r["term_id"]] = (r["term_label"], split(r["parent_ids"]),
                                                split(r["ancestor_ids"]))
    return terms


def species_map(datasets):
    """{organism id: (species id, species label)}: subspecies folded into their species."""
    by_label = {}
    for d in datasets:
        for o in d["organism"]:
            by_label.setdefault(o["label"], o["ontology_term_id"])
    fold, missing = {}, []
    for label, oid in by_label.items():
        words = label.split()
        species = " ".join(words[:2])
        if len(words) > 2 and species in by_label:
            fold[oid] = (by_label[species], species)
        else:
            if len(words) > 2:
                missing.append(label)
            fold[oid] = (oid, label)
    if missing:
        print("warning: subspecies kept, their species is in no dataset: " + ", ".join(missing))
    return fold


def term_values(dataset, field):
    """[(id, label)] of an ontology-term field; unknown/na/empty -> 'na'."""
    out = {}
    for o in dataset.get(field) or []:
        tid = o.get("ontology_term_id") or f"label:{o.get('label')}"
        if tid in UNKNOWN_IDS:
            tid = NA
        out.setdefault(tid, NA_LABEL if tid == NA else o.get("label") or tid)
    return list(out.items()) or [(NA, NA_LABEL)]


def title_of(d):
    if d.get("title"):
        return d["title"]
    if d.get("collection_name"):
        return d["collection_name"]
    ref = d["source_reference_id"]
    return f"{d['collection_doi_label']} ({ref})" if d.get("collection_doi_label") \
        else f"{d['source']['name']} {ref}"


def publications(datasets):
    """Publication (or study) of each dataset; see 'publication' in the module docstring.

    Returns ({dataset_id: publication id}, {publication id: (label, url)}, {ignored collection:
    its DOIs})."""
    def collection(d):
        c, slug = d.get("collection_id"), d["source"]["slug"]
        return (slug, c) if c and slug in COLLECTION_RESOURCES and c != d["source_reference_id"] else None

    coll_dois = collections.defaultdict(set)
    for d in datasets:
        if collection(d):
            coll_dois[collection(d)].add(norm_doi(d.get("collection_doi")))
    ignored = {c: sorted(x for x in dois if x) for c, dois in coll_dois.items()
               if len(dois - {None}) > 1}

    pub_of, members = {}, collections.defaultdict(list)
    for d in datasets:
        c = collection(d)
        if c in ignored:
            c = None
        doi = norm_doi(d.get("collection_doi"))
        if not doi and c and coll_dois[c] - {None}:
            doi = next(iter(coll_dois[c] - {None}))  # DOI of the other datasets of the collection
        key = f"doi:{doi}" if doi else f"collection:{c[0]}:{c[1]}" if c else f"dataset:{d['dataset_id']}"
        pub_of[d["dataset_id"]] = key
        members[key].append(d)

    info = {}
    for key, ds in members.items():
        def most(field):
            values = collections.Counter(d.get(field) for d in ds if d.get(field))
            return values.most_common(1)[0][0] if values else None
        if key.startswith("doi:"):
            doi = clean_doi(most("collection_doi")) or key[4:]
            info[key] = (most("collection_doi_label") or most("collection_name") or doi,
                         f"https://doi.org/{doi}")
        elif key.startswith("collection:"):
            cid = ds[0]["collection_id"]
            url = next((d["source_url"] for d in ds if cid in (d.get("source_url") or "")), "")
            info[key] = (most("collection_name") or f"{ds[0]['source']['name']} collection {cid}", url)
        else:
            info[key] = (title_of(ds[0]), "")
    by_label = collections.defaultdict(list)
    for key, (label, _) in info.items():
        by_label[label].append(key)
    for label, keys in by_label.items():
        if len(keys) > 1:  # different publications with the same label
            for key in keys:
                info[key] = (f"{label} [{key.split(':', 1)[1]}]", info[key][1])
    return pub_of, info, ignored


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--datasets", type=Path, required=True,
                    help="merged datasets JSON written by merge_datasets.py")
    ap.add_argument("--remappings", type=Path, required=True,
                    help="term remapping TSV written by make_term_remappings.py")
    ap.add_argument("--ancestors", type=Path, required=True,
                    help="term ancestor TSV written by make_term_ancestors.py")
    ap.add_argument("--subsets", type=Path, required=True,
                    help="term subset TSV written by make_term_subsets.py")
    ap.add_argument("--out-datasets", type=Path, required=True, help="dataset Parquet to write")
    ap.add_argument("--out-terms", type=Path, required=True, help="term Parquet to write")
    args = ap.parse_args()

    datasets = json.load(args.datasets.open(encoding="utf-8"))
    remap = read_remappings(args.remappings)
    hierarchy = read_ancestors(args.ancestors)
    fold = species_map(datasets)
    pub_of, pub_info, ignored_collections = publications(datasets)
    with args.subsets.open(encoding="utf-8", newline="") as f:
        subsets = {(r["factor"], r["term_id"]): r["subsets"].split(";")
                   for r in csv.DictReader(f, delimiter="\t")}

    labels = collections.defaultdict(dict)  # field -> {value id: label}
    missing = collections.Counter()
    rows = []
    for d in datasets:
        row = {"id": d["dataset_id"], "title": title_of(d)}
        species = {}
        for o in d["organism"]:
            sid, slabel = fold[o["ontology_term_id"]]
            species.setdefault(sid, slabel)
        row["organism"] = list(species)
        labels["organism"].update(species)

        for factor in FACTORS:
            # source labels, for terms missing from the ancestor TSV (unresolved remappings)
            annotated = {o["ontology_term_id"]: o["label"] for o in d.get(factor) or []}
            terms = remapped_terms(d, factor, remap)
            for t in terms:
                if t not in hierarchy[factor]:
                    missing[f"{factor}:{t}"] += 1
                    labels[factor].setdefault(t, annotated.get(t, t))
            row[factor] = terms

        for field in TERM_FIELDS:
            values = term_values(d, field)
            row[field] = [v for v, _ in values]
            for v, label in values:
                labels[field].setdefault(v, label)

        row["suspension_type"] = list(dict.fromkeys(d.get("suspension_type") or [NA]))
        types = list(dict.fromkeys(o["tissue_type"] for o in d.get("tissue") or []
                                   if o.get("tissue_type")))
        row["tissue_type"] = types or [NA]
        for field in ("suspension_type", "tissue_type"):
            for v in row[field]:
                labels[field].setdefault(v, NA_LABEL if v == NA else v)

        src = {}
        for s in d["sources"]:
            src.setdefault(s["source"]["slug"], (s["source"]["name"], s.get("explorer_url") or ""))
        primary = d["source"]["slug"]  # primary resource first, its URL is the dataset's
        src[primary] = (src[primary][0], d.get("explorer_url") or src[primary][1])
        order = [primary] + [s for s in src if s != primary]
        row["source"] = order
        row["source_explorer_urls"] = [src[s][1] for s in order]
        labels["source"].update({s: src[s][0] for s in order})

        pub = pub_of[d["dataset_id"]]
        row["publication"] = [pub]
        row["publication_url"] = pub_info[pub][1]
        labels["publication"][pub] = pub_info[pub][0]

        row["flags"] = [f for f, on in (("merged", d.get("merged")),
                                        ("multi_species", len(species) > 1)) if on]
        labels["flags"].update({f: FLAG_LABELS[f] for f in row["flags"]})

        row["cell_count"] = int(d.get("cell_count") or 0)
        row["explorer_url"] = d.get("explorer_url") or ""
        row["h5ad_url"] = next((a["url"] for a in d.get("assets") or []
                                if str(a.get("filetype", "")).upper() == "H5AD"), "")
        row["doi"] = d.get("collection_doi") or ""
        row["doi_label"] = d.get("collection_doi_label") or ""
        row["collection_name"] = d.get("collection_name") or ""
        rows.append(row)

    # Term table: the hierarchy of the factors, plain values for the other fields.
    term_rows = []
    for factor in FACTORS:
        for tid, (label, parents, ancestors) in hierarchy[factor].items():
            term_rows.append({"field": factor, "id": tid, "label": label,
                              "parent_ids": parents, "ancestor_ids": ancestors,
                              "subsets": subsets.get((factor, tid), [])})
    for field, values in labels.items():
        for tid, label in values.items():
            if field in FACTORS and tid in hierarchy[field]:
                continue
            term_rows.append({"field": field, "id": tid, "label": label,
                              "parent_ids": [], "ancestor_ids": [], "subsets": []})

    # Text search: title, ids and the labels of the values (not of their ancestors).
    term_label = {(r["field"], r["id"]): r["label"] for r in term_rows}
    for row in rows:
        bits = [row["title"], row["id"], row["doi"], row["doi_label"], row["collection_name"]]
        for field in ("organism", *FACTORS, *TERM_FIELDS, "source", "publication"):
            bits += [term_label[(field, v)] for v in row[field] if v != NA]
        row["search_blob"] = " | ".join(b for b in bits if b)

    args.out_datasets.parent.mkdir(parents=True, exist_ok=True)
    args.out_terms.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    table = pa.table({c: pa.array([r[c] for r in rows],
                                  type=pa.list_(pa.string()) if c in LIST_COLUMNS
                                  else pa.int64() if c == "cell_count" else pa.string())
                      for c in columns})
    pq.write_table(table, args.out_datasets, compression="zstd")
    terms = pa.table({c: pa.array([r[c] for r in term_rows],
                                  type=pa.list_(pa.string()) if c.endswith("_ids")
                                  or c == "subsets" else pa.string())
                      for c in ("field", "id", "label", "parent_ids", "ancestor_ids", "subsets")})
    pq.write_table(terms, args.out_terms, compression="zstd")

    if missing:
        print(f"warning: {len(missing)} annotated terms not in the ancestor TSV (no roll-up): "
              + ", ".join(sorted(missing)))
    print(f"{len(rows)} datasets -> {args.out_datasets} "
          f"({args.out_datasets.stat().st_size / 1e6:.2f} MB)")
    print(f"{len(term_rows)} terms -> {args.out_terms} "
          f"({args.out_terms.stat().st_size / 1e6:.2f} MB)")
    hidden = collections.Counter(r["field"] for r in term_rows if r["subsets"])
    for field, n in collections.Counter(r["field"] for r in term_rows).items():
        print(f"  {field:<18} {n:>5} values" + (f", {hidden[field]} in subsets" if hidden[field] else ""))
    kinds = collections.Counter(p.split(":", 1)[0] for p in pub_info)
    in_kind = collections.Counter(p.split(":", 1)[0] for p in pub_of.values())
    print(f"  publications: {len(pub_info)} for {len(rows)} datasets: "
          + ", ".join(f"{k} {kinds[k]} ({in_kind[k]} datasets)" for k in ("doi", "collection", "dataset")))
    for (slug, cid), dois in sorted(ignored_collections.items()):
        print(f"  collection {slug}:{cid} ignored for publications: its datasets have {len(dois)} DOIs")
    unknown = set(subsets) - {(r["field"], r["id"]) for r in term_rows}
    if unknown:
        print(f"warning: {len(unknown)} terms of the subset TSV are not in the ancestor TSV")


if __name__ == "__main__":
    main()
