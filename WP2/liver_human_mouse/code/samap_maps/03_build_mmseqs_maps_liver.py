#!/usr/bin/env python3
"""Continue 02_build_samap_fastas.py: build mouse/zebrafish MMseqs2 maps.

Reuses the existing selected cDNA sequences and SAMap gene names unchanged.
Run in your SAMap environment after installing MMseqs2:
    conda install -c conda-forge -c bioconda mmseqs2
    python scripts/03_build_mmseqs_maps.py

Documentation: https://github.com/soedinglab/MMseqs2/wiki#translated-sequence-searching

Identifier fix: export full FASTA headers, then take their first token. MMseqs2
accession parsing can otherwise truncate custom names such as 'ucp3' to 'p3'.
On failure, working files are retained in ROOT/tmp instead of being deleted.
"""

from pathlib import Path
import csv
import math
import shutil
import subprocess
import tempfile


# Same project directory and FASTA files as your scripts/notebook.
# LIVER ADAPTATION 2026-11-11: copied from fbastian/BH2026-23-single-cell WP2/scripts (UPSTREAM_COMMIT.txt);
# changes: ROOT, species pair hs/mm and file names, THREADS. Search parameters unchanged.
ROOT = Path("/mnt/masld_vault/cross_species/samap_maps/liver_hsmm")
FASTAS = {
    "hs": ROOT / "refs/samap/human.GRCh38.GENCODEv48.samap.fa",
    "mm": ROOT / "refs/samap/mouse.GRCm38.GENCODEvM18.samap.fa",
}
OUT = ROOT / "maps_mmseqs" / "hsmm"
THREADS = 20  # TotalVI shares the VM
EVALUE = 1e-6
SENSITIVITY = 7.5
# Prefilter candidates per translated query; increase for a sensitivity check.
MAX_SEQS = 10000
# 'query,target' invoke MMseqs2 accession parsing; full headers bypass that parsing.
FIELDS = "qheader,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"


def run(*arguments):
    subprocess.run([str(arg) for arg in arguments], check=True)


def fasta_names(path):
    """Read the first header token, exactly as used by your existing BLAST maps."""
    names = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    if not names or len(names) != len(set(names)):
        raise ValueError(f"Empty FASTA or duplicate gene names: {path}")
    return set(names)


def keep_best_alignment(raw_path, output_path, query_names, target_names):
    """Collapse ORF/frame hits to the best alignment per directed gene pair.

    Keep multiple target genes per query. SAMap handles reciprocity itself.
    Preserve the selected alignment's real E-value and bit score.
    """
    best = {}
    with raw_path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 12:
                raise ValueError("Expected 12 alignment columns")
            # Script 02 writes: >samap_gene gene_id=... transcript_id=... length=...
            # Keep the exact first token, including colons and '__EnsemblGeneID'.
            if not row[0].split() or not row[1].split():
                raise ValueError("Empty query or target FASTA header in alignment output")
            query, target = row[0].split()[0], row[1].split()[0]
            if query not in query_names or target not in target_names:
                raise ValueError(f"Full FASTA headers do not match the input gene names: {query}, {target}")
            row[0], row[1] = query, target
            evalue, bits = float(row[10]), float(row[11])
            if not math.isfinite(evalue) or not math.isfinite(bits) or evalue < 0 or bits <= 0:
                raise ValueError("Invalid E-value or bit score")
            if evalue > EVALUE:
                continue
            pair = (query, target)
            # Highest bit score, then lowest E-value; deterministic exact ties.
            candidate = (bits, -evalue, tuple(row))
            if pair not in best or candidate > best[pair]:
                best[pair] = candidate
    if not best:
        raise ValueError(f"No significant gene pairs in {raw_path}")
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(best[pair][2] for pair in sorted(best))
    print(f"{output_path.name}: {len(best):,} gene pairs", flush=True)


def main():
    if shutil.which("mmseqs") is None:
        raise RuntimeError("MMseqs2 is missing; install it with the command in this script's header")
    if OUT.exists():
        raise FileExistsError(f"Output already exists: {OUT}; choose a different OUT for a new run")
    names = {sid: fasta_names(path) for sid, path in FASTAS.items()}
    temporary_root = ROOT / "tmp"
    temporary_root.mkdir(parents=True, exist_ok=True)

    work = Path(tempfile.mkdtemp(prefix="mmseqs_hs_mm_", dir=temporary_root))
    print(f"Working directory: {work}", flush=True)
    try:
        for sid, path in FASTAS.items():
            # Explicit nucleotide databases: do not interpret cDNA as protein.
            run("mmseqs", "createdb", path, work / sid, "--dbtype", 2)

        for query, target in [("hs", "mm"), ("mm", "hs")]:
            result = work / f"{query}_{target}_result"
            raw = work / f"{query}_{target}_raw.tsv"
            run("mmseqs", "search", work / query, work / target, result,
                work / f"tmp_{query}_{target}",
                "--search-type", 2, "-s", SENSITIVITY, "-e", EVALUE,
                "--min-length", 30, "--orf-start-mode", 1, "--translation-table", 1,
                "--max-seqs", MAX_SEQS, "--max-accept", 1000000,
                "--max-rejected", 1000000, "-a", "--threads", THREADS)
            run("mmseqs", "convertalis", work / query, work / target, result, raw,
                "--search-type", 2, "--format-output", FIELDS, "--threads", THREADS)
            keep_best_alignment(raw, work / f"{query}_to_{target}.txt",
                                names[query], names[target])

        # Publish the two finished maps together, separately from your BLAST maps.
        OUT.mkdir(parents=True)
        for filename in ("hs_to_mm.txt", "mm_to_hs.txt"):
            shutil.copyfile(work / filename, OUT / filename)
        with (OUT / "mmseqs_version.txt").open("w") as handle:
            subprocess.run(["mmseqs", "version"], check=True, stdout=handle)
    except BaseException:
        print(f"Run interrupted or failed. Working files retained in: {work}", flush=True)
        raise
    else:
        shutil.rmtree(work)
    print(f"Finished. SAMap f_maps: {OUT.parent}/")


if __name__ == "__main__":
    main()
