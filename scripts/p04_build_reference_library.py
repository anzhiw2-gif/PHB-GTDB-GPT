"""Build and validate the P04 reference-library manifest.

This stage keeps the biological evidence boundary explicit:
reference seeds, family labels, accession provenance, and retrieval metadata
are tracked as auditable rows rather than inferred from downstream hits.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from datetime import date
from pathlib import Path


CANONICAL_EVIDENCE_LEVELS = ("E1", "E2", "E3", "Excluded")
CANONICAL_TAXONOMIC_DOMAINS = ("Bacteria", "Archaea")
CANONICAL_REFERENCE_LIBRARIES = (
    "bacteria_high_confidence",
    "archaea_literature_supported",
)
CANONICAL_FAMILY_CATEGORIES = (
    "intracellular_phaZ_no_lipase_box",
    "extracellular_scl_pha_dep_type_I",
    "extracellular_scl_pha_dep_type_II",
    "phaZ7_like",
    "phaZd_like",
    "rhodospirillum_periplasmic_like",
    "intracellular_mcl_pha_dep",
    "extracellular_mcl_pha_dep",
    "tigr02240_aromatic_pha_related",
    "archaeal_patatin_like_pha_dep",
    "auxiliary_mobilization_context",
)
CANONICAL_SOURCE_DATABASES = (
    "UniProtKB",
    "NCBI Protein",
    "NCBI Nucleotide",
    "PDB",
    "Legacy verified",
)
CANONICAL_SEQUENCE_FORMATS = ("fasta", "faa")

REQUIRED_FIELDS = (
    "seed_id",
    "reference_library",
    "taxonomic_domain",
    "family_category",
    "seed_name",
    "evidence_level",
    "source_database",
    "source_accession",
    "organism",
    "taxon_id",
    "retrieval_date",
    "sequence_format",
    "sequence_length_aa",
    "sequence_path",
)

OPTIONAL_FIELDS = (
    "family_label",
    "source_release",
    "source_version",
    "source_url",
    "retrieval_method",
    "retrieval_query",
    "retrieval_endpoint",
    "retrieval_batch_id",
    "retrieval_log_path",
    "accession_version",
    "taxon",
    "database_version",
    "doi",
    "pmid",
    "pmcid",
    "reference_title",
    "reference_year",
    "exclusion_reason",
    "record_kind",
    "literature_support_scope",
    "supporting_sources",
    "supporting_accessions",
    "supporting_notes",
    "notes",
)

OUTPUT_FILENAME = "reference_library.normalized.tsv"
BACTERIA_OUTPUT_FILENAME = "reference_library.bacteria.normalized.tsv"
ARCHAEA_OUTPUT_FILENAME = "reference_library.archaea.normalized.tsv"
SUMMARY_FILENAME = "reference_library_summary.tsv"


def load_reference_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing_columns = [field for field in REQUIRED_FIELDS if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            validate_reference_row(normalized, path=path, line_number=line_number)
            rows.append(normalized)
    return rows


def validate_reference_row(row: dict[str, str], path: Path, line_number: int) -> None:
    for field in REQUIRED_FIELDS:
        if not row.get(field):
            raise ValueError(f"{path}:{line_number} missing required value for {field}")

    if row["family_category"] not in CANONICAL_FAMILY_CATEGORIES:
        allowed = ", ".join(CANONICAL_FAMILY_CATEGORIES)
        raise ValueError(
            f"{path}:{line_number} invalid family_category {row['family_category']!r}; allowed: {allowed}"
        )

    if row["reference_library"] not in CANONICAL_REFERENCE_LIBRARIES:
        allowed = ", ".join(CANONICAL_REFERENCE_LIBRARIES)
        raise ValueError(
            f"{path}:{line_number} invalid reference_library {row['reference_library']!r}; allowed: {allowed}"
        )

    if row["taxonomic_domain"] not in CANONICAL_TAXONOMIC_DOMAINS:
        allowed = ", ".join(CANONICAL_TAXONOMIC_DOMAINS)
        raise ValueError(
            f"{path}:{line_number} invalid taxonomic_domain {row['taxonomic_domain']!r}; allowed: {allowed}"
        )

    if row["evidence_level"] not in CANONICAL_EVIDENCE_LEVELS:
        allowed = ", ".join(CANONICAL_EVIDENCE_LEVELS)
        raise ValueError(
            f"{path}:{line_number} invalid evidence_level {row['evidence_level']!r}; allowed: {allowed}"
        )

    validate_library_boundary(row, path=path, line_number=line_number)

    if row["source_database"] not in CANONICAL_SOURCE_DATABASES:
        allowed = ", ".join(CANONICAL_SOURCE_DATABASES)
        raise ValueError(
            f"{path}:{line_number} invalid source_database {row['source_database']!r}; allowed: {allowed}"
        )

    try:
        date.fromisoformat(row["retrieval_date"])
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} invalid retrieval_date {row['retrieval_date']!r}; expected YYYY-MM-DD"
        ) from error

    if not row["taxon_id"].isdigit():
        raise ValueError(f"{path}:{line_number} taxon_id must be a numeric NCBI taxon ID")

    sequence_length = _parse_sequence_length(row["sequence_length_aa"], path=path, line_number=line_number)
    if sequence_length < 1:
        raise ValueError(f"{path}:{line_number} sequence_length_aa must be a positive integer")

    if row["sequence_format"] not in CANONICAL_SEQUENCE_FORMATS:
        allowed = ", ".join(CANONICAL_SEQUENCE_FORMATS)
        raise ValueError(
            f"{path}:{line_number} invalid sequence_format {row['sequence_format']!r}; allowed: {allowed}"
        )

    if row["evidence_level"] == "Excluded" and not row.get("exclusion_reason"):
        raise ValueError(f"{path}:{line_number} exclusion_reason is required when evidence_level is Excluded")


def validate_library_boundary(row: dict[str, str], path: Path, line_number: int) -> None:
    expected_library = {
        "Bacteria": "bacteria_high_confidence",
        "Archaea": "archaea_literature_supported",
    }[row["taxonomic_domain"]]
    if row["reference_library"] != expected_library:
        raise ValueError(
            f"{path}:{line_number} taxonomic_domain {row['taxonomic_domain']!r} "
            f"requires reference_library {expected_library!r}"
        )

    if row["evidence_level"] == "Excluded":
        return

    if row["reference_library"] == "bacteria_high_confidence" and row["evidence_level"] not in {"E1", "E2"}:
        raise ValueError(
            f"{path}:{line_number} bacteria_high_confidence rows require E1 or E2 evidence"
        )

    if row["reference_library"] == "archaea_literature_supported":
        if row["evidence_level"] not in {"E1", "E2"}:
            raise ValueError(
                f"{path}:{line_number} archaea_literature_supported rows require E1 or E2 evidence"
            )
        if not (row.get("pmid") or row.get("doi") or row.get("pmcid")):
            raise ValueError(
                f"{path}:{line_number} archaea_literature_supported rows require PMID, DOI, or PMCID"
            )
        if not row.get("literature_support_scope"):
            raise ValueError(
                f"{path}:{line_number} literature_support_scope is required for archaeal literature-supported rows"
            )


def normalize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            row["reference_library"],
            row["family_category"],
            _evidence_rank(row["evidence_level"]),
            row["seed_id"],
        ),
    )


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    by_library = Counter(row["reference_library"] for row in rows)
    by_domain = Counter(row["taxonomic_domain"] for row in rows)
    by_family = Counter(row["family_category"] for row in rows)
    by_evidence = Counter(row["evidence_level"] for row in rows)
    summary: list[dict[str, object]] = []

    for reference_library in CANONICAL_REFERENCE_LIBRARIES:
        summary.append(
            {
                "kind": "reference_library",
                "name": reference_library,
                "count": by_library.get(reference_library, 0),
            }
        )

    for taxonomic_domain in CANONICAL_TAXONOMIC_DOMAINS:
        summary.append(
            {
                "kind": "taxonomic_domain",
                "name": taxonomic_domain,
                "count": by_domain.get(taxonomic_domain, 0),
            }
        )

    for family_category in CANONICAL_FAMILY_CATEGORIES:
        summary.append(
            {
                "kind": "family",
                "name": family_category,
                "count": by_family.get(family_category, 0),
            }
        )

    for evidence_level in CANONICAL_EVIDENCE_LEVELS:
        summary.append(
            {
                "kind": "evidence",
                "name": evidence_level,
                "count": by_evidence.get(evidence_level, 0),
            }
        )

    return summary


def build_reference_library(manifest_path: Path, outdir: Path) -> dict[str, Path]:
    rows = load_reference_manifest(manifest_path)
    normalized = normalize_rows(rows)
    outdir.mkdir(parents=True, exist_ok=True)

    normalized_path = outdir / OUTPUT_FILENAME
    bacteria_path = outdir / BACTERIA_OUTPUT_FILENAME
    archaea_path = outdir / ARCHAEA_OUTPUT_FILENAME
    summary_path = outdir / SUMMARY_FILENAME

    write_tsv(normalized_path, normalized, REQUIRED_FIELDS + OPTIONAL_FIELDS)
    write_tsv(
        bacteria_path,
        [row for row in normalized if row["reference_library"] == "bacteria_high_confidence"],
        REQUIRED_FIELDS + OPTIONAL_FIELDS,
    )
    write_tsv(
        archaea_path,
        [row for row in normalized if row["reference_library"] == "archaea_literature_supported"],
        REQUIRED_FIELDS + OPTIONAL_FIELDS,
    )
    write_tsv(summary_path, summarize_rows(rows), ("kind", "name", "count"))

    return {
        "normalized": normalized_path,
        "bacteria": bacteria_path,
        "archaea": archaea_path,
        "summary": summary_path,
    }


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build the P04 reference-library manifest.")
    parser.add_argument("--manifest", required=True, type=Path, help="Input curated seed manifest TSV")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for normalized assets")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    outputs = build_reference_library(args.manifest, args.outdir)
    print(f"Normalized manifest written: {outputs['normalized']}")
    print(f"Summary written: {outputs['summary']}")
    return 0


def _evidence_rank(level: str) -> int:
    return CANONICAL_EVIDENCE_LEVELS.index(level)


def _parse_sequence_length(value: str, path: Path, line_number: int) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"{path}:{line_number} sequence_length_aa must be an integer") from error


if __name__ == "__main__":
    raise SystemExit(main())
