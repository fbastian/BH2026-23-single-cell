from pathlib import Path
import gzip
import re
import textwrap

import pandas as pd
import scanpy as sc


# ============================================================
# CONFIGURATION
# ============================================================

# LIVER ADAPTATION (2026-11-11): copied from fbastian/BH2026-23-single-cell WP2/scripts (commit in
# UPSTREAM_COMMIT.txt). Changes: paths, species = human HLiCA (GENCODE v48) + mouse Guilliams 2022 (GENCODE vM18 =
# Ensembl 93), gene fields of our objects, GENCODE parser accepts human IDs. Naming + transcript choice unchanged.
ROOT = Path("/mnt/masld_vault/cross_species/samap_maps/liver_hsmm")

RAW = ROOT / "refs" / "raw"
OUT = ROOT / "refs" / "samap"
MAPS = ROOT / "refs" / "mappings"

OUT.mkdir(parents=True, exist_ok=True)
MAPS.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# MODIFY ONLY THESE PATH
# ------------------------------------------------------------

HUMAN_H5AD = Path(
    "/mnt/masld_vault/geo/cxg-059202e1-1f1b-483f-9151-f3a25a380c39/HLiCA_human_liver_all_cells.h5ad"
)

MOUSE_H5AD = Path(
    "/mnt/masld_vault/cross_species/guilliams2022/01_raw/guilliams2022_fig1b_raw.h5ad"
)


CONFIG = {
    "hs": {
        "name": "human",
        "h5ad": HUMAN_H5AD,

        "fasta": (
            RAW
            / "gencode.v48.transcripts.fa.gz"
        ),

        "fasta_format": "gencode",

        "output": (
            OUT
            / "human.GRCh38.GENCODEv48.samap.fa"
        ),
    },

    "mm": {
        "name": "mouse",
        "h5ad": MOUSE_H5AD,

        "fasta": (
            RAW
            / "gencode.vM18.transcripts.fa.gz"
        ),

        "fasta_format": "gencode",

        "output": (
            OUT
            / "mouse.GRCm38.GENCODEvM18.samap.fa"
        ),
    },
}


# ============================================================
# UTILITIES
# ============================================================

def strip_version(identifier):
    """
    ENSDARG00000012345.7
          ->
    ENSDARG00000012345
    """

    identifier = str(identifier)

    return re.sub(
        r"\.\d+$",
        "",
        identifier,
    )


def fasta_reader(path):
    """
    Streaming FASTA reader for .fa.gz.
    """

    header = None
    sequence = []

    with gzip.open(path, "rt") as handle:

        for line in handle:

            line = line.rstrip()

            if line.startswith(">"):

                if header is not None:
                    yield header, "".join(sequence)

                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        if header is not None:
            yield header, "".join(sequence)


# ============================================================
# EXTRACT TRANSCRIPT / GENE IDS
# ============================================================

def parse_ensembl_header(header):
    """
    Expected Ensembl style:

    >ENSDART... chromosome:... gene:ENSDARG...
    """

    transcript_id = header.split()[0]

    match = re.search(
        r"(?:^|\s)gene:([^\s]+)",
        header,
    )

    if match is None:
        return None, None

    gene_id = strip_version(
        match.group(1)
    )

    transcript_id = strip_version(
        transcript_id
    )

    return transcript_id, gene_id


def parse_gencode_header(header):
    """
    GENCODE transcript FASTA style:

    >ENSMUST....|ENSMUSG....|...
    """

    first = header.split()[0]

    fields = first.split("|")

    if len(fields) < 2:
        return None, None

    transcript_id = strip_version(
        fields[0]
    )

    gene_id = strip_version(
        fields[1]
    )

    if not gene_id.startswith(
        ("ENSMUSG", "ENSG")
    ):
        return None, None

    return transcript_id, gene_id


# ============================================================
# CREATE DATASET GENE TABLE
# ============================================================

def load_gene_table(species, h5ad_file):

    print()
    print("=" * 70)
    print("Reading H5AD:", species)
    print("=" * 70)

    # backed mode:
    # we only need .var, not the huge expression matrix
    adata = sc.read_h5ad(
        h5ad_file,
        backed="r",
    )

    # liver objects: var_names = Ensembl IDs; symbols in var['feature_name'] (HLiCA, CxG) / var['symbol'] (mouse)
    if species == "hs":

        gene_id = pd.Series(
            adata.var_names.astype(str)
        )

        gene_symbol = pd.Series(
            adata.var["feature_name"].astype(str).values
        )

    elif species == "mm":

        gene_id = pd.Series(
            adata.var_names.astype(str)
        )

        gene_symbol = pd.Series(
            adata.var[
                "symbol"
            ].astype(str).values
        )

    else:
        raise ValueError(species)

    adata.file.close()

    gene_id = gene_id.map(
        strip_version
    )

    gene_symbol = (
        gene_symbol
        .astype(str)
        .str.strip()
    )


    # --------------------------------------------------------
    # Sanity checks
    # --------------------------------------------------------

    if gene_id.isna().any():
        raise RuntimeError(
            f"{species}: missing gene IDs"
        )

    if gene_symbol.isna().any():
        raise RuntimeError(
            f"{species}: missing gene symbols"
        )

    if (gene_symbol == "").any():
        raise RuntimeError(
            f"{species}: empty gene symbols"
        )


    # FASTA IDs are terminated by whitespace.
    # Gene symbols should therefore not contain spaces.
    has_space = gene_symbol.str.contains(
        r"\s",
        regex=True,
    )

    if has_space.any():

        print(
            gene_symbol[
                has_space
            ].head(20)
        )

        raise RuntimeError(
            f"{species}: gene symbols "
            "containing whitespace found."
        )


    # --------------------------------------------------------
    # Detect ambiguous symbols
    # --------------------------------------------------------

    duplicated_symbol = (
        gene_symbol
        .duplicated(
            keep=False
        )
    )


    samap_gene = gene_symbol.copy()


    # IMPORTANT:
    # every member of an ambiguous symbol group
    # gets its Ensembl Gene ID.
    samap_gene.loc[
        duplicated_symbol
    ] = (
        gene_symbol.loc[
            duplicated_symbol
        ]
        + "__"
        + gene_id.loc[
            duplicated_symbol
        ]
    )


    if not pd.Index(
        samap_gene
    ).is_unique:

        dup = samap_gene[
            samap_gene.duplicated(
                keep=False
            )
        ]

        raise RuntimeError(
            "SAMap names remain duplicated:\n"
            + str(dup.head(30))
        )


    table = pd.DataFrame(
        {
            "ensembl_gene_id":
                gene_id.values,

            "gene_symbol":
                gene_symbol.values,

            "samap_gene":
                samap_gene.values,

            "symbol_ambiguous":
                duplicated_symbol.values,
        }
    )


    table.to_csv(
        MAPS
        / f"{species}.gene_names.tsv",
        sep="\t",
        index=False,
    )


    print(
        "Genes:",
        len(table)
    )

    print(
        "Ambiguous genes:",
        table[
            "symbol_ambiguous"
        ].sum()
    )

    print(
        "Final SAMap names unique:",
        table[
            "samap_gene"
        ].is_unique
    )


    if table[
        "symbol_ambiguous"
    ].any():

        print()
        print(
            "Examples of ambiguous symbols:"
        )

        print(
            table.loc[
                table[
                    "symbol_ambiguous"
                ],
                [
                    "ensembl_gene_id",
                    "gene_symbol",
                    "samap_gene",
                ]
            ].head(20)
        )


    return table


# ============================================================
# BUILD ONE-SEQUENCE-PER-GENE FASTA
# ============================================================

def build_fasta(
    species,
    gene_table,
    input_fasta,
    output_fasta,
    fasta_format,
):

    print()
    print("=" * 70)
    print(
        "Building FASTA:",
        species
    )
    print("=" * 70)


    gene_to_name = dict(
        zip(
            gene_table[
                "ensembl_gene_id"
            ],
            gene_table[
                "samap_gene"
            ],
        )
    )


    if fasta_format == "ensembl":
        parser = parse_ensembl_header

    elif fasta_format == "gencode":
        parser = parse_gencode_header

    else:
        raise ValueError(
            fasta_format
        )


    # --------------------------------------------------------
    # best[gene_id] =
    # {
    #    transcript_id,
    #    sequence
    # }
    #
    # Criterion:
    # longest cDNA transcript
    # --------------------------------------------------------

    best = {}

    total_records = 0
    records_with_gene_id = 0
    records_matching_dataset = 0


    for header, sequence in fasta_reader(
        input_fasta
    ):

        total_records += 1

        transcript_id, gene_id = parser(
            header
        )

        if gene_id is None:
            continue

        records_with_gene_id += 1


        # Gene not used by the scRNA-seq dataset
        if gene_id not in gene_to_name:
            continue

        records_matching_dataset += 1


        current = best.get(
            gene_id
        )


        replace = False

        if current is None:
            replace = True

        elif len(sequence) > len(
            current["sequence"]
        ):
            replace = True

        elif (
            len(sequence)
            == len(
                current["sequence"]
            )
            and transcript_id
            < current["transcript_id"]
        ):
            # deterministic tie breaking
            replace = True


        if replace:

            best[gene_id] = {
                "transcript_id":
                    transcript_id,

                "sequence":
                    sequence,
            }


    print(
        "Total FASTA records:",
        total_records
    )

    print(
        "Records with parsed gene ID:",
        records_with_gene_id
    )

    print(
        "Records belonging to dataset genes:",
        records_matching_dataset
    )


    # --------------------------------------------------------
    # Write output in DATASET order
    # --------------------------------------------------------

    audit_rows = []

    n_written = 0


    with open(
        output_fasta,
        "w"
    ) as out:

        for row in gene_table.itertuples(
            index=False
        ):

            gene_id = (
                row.ensembl_gene_id
            )

            gene_symbol = (
                row.gene_symbol
            )

            samap_gene = (
                row.samap_gene
            )


            selected = best.get(
                gene_id
            )


            if selected is None:

                audit_rows.append(
                    {
                        "ensembl_gene_id":
                            gene_id,

                        "gene_symbol":
                            gene_symbol,

                        "samap_gene":
                            samap_gene,

                        "transcript_id":
                            None,

                        "length":
                            None,

                        "found":
                            False,
                    }
                )

                continue


            transcript_id = selected[
                "transcript_id"
            ]

            sequence = selected[
                "sequence"
            ]


            # ------------------------------------------------
            # IMPORTANT:
            #
            # first FASTA token is the exact SAMap gene name.
            #
            # Additional metadata comes AFTER whitespace,
            # so BLAST qseqid remains samap_gene.
            # ------------------------------------------------

            out.write(
                f">{samap_gene} "
                f"gene_id={gene_id} "
                f"transcript_id={transcript_id} "
                f"length={len(sequence)}\n"
            )


            for line in textwrap.wrap(
                sequence,
                width=80,
            ):
                out.write(
                    line + "\n"
                )


            n_written += 1


            audit_rows.append(
                {
                    "ensembl_gene_id":
                        gene_id,

                    "gene_symbol":
                        gene_symbol,

                    "samap_gene":
                        samap_gene,

                    "transcript_id":
                        transcript_id,

                    "length":
                        len(sequence),

                    "found":
                        True,
                }
            )


    audit = pd.DataFrame(
        audit_rows
    )


    audit_file = (
        MAPS
        / f"{species}.selected_transcripts.tsv"
    )

    audit.to_csv(
        audit_file,
        sep="\t",
        index=False,
    )


    missing = audit.loc[
        ~audit["found"]
    ]


    missing.to_csv(
        MAPS
        / f"{species}.missing_from_transcriptome.tsv",
        sep="\t",
        index=False,
    )


    print()
    print(
        "Genes written:",
        n_written
    )

    print(
        "Dataset genes:",
        len(gene_table)
    )

    print(
        "Genes missing from FASTA:",
        len(missing)
    )


    coverage = (
        100
        * n_written
        / len(gene_table)
    )


    print(
        f"Coverage: {coverage:.2f}%"
    )

    print()
    print(
        "FASTA:",
        output_fasta
    )

    print(
        "Audit table:",
        audit_file
    )


# ============================================================
# MAIN
# ============================================================

for species, cfg in CONFIG.items():

    gene_table = load_gene_table(
        species,
        cfg["h5ad"],
    )


    build_fasta(
        species=species,

        gene_table=gene_table,

        input_fasta=cfg["fasta"],

        output_fasta=cfg["output"],

        fasta_format=cfg[
            "fasta_format"
        ],
    )


print()
print("=" * 70)
print("DONE")
print("=" * 70)
