"""Prepare deterministic MAFFT and HMMER commands for P05 family profiles.

This step only materializes a command manifest. It does not execute MAFFT,
hmmbuild, or downstream calibration. Sequence evidence remains a family-model
input and is not a phenotype claim.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

try:  # pragma: no cover - direct script execution from repo root
    from scripts.p05_plan_family_profiles import (
        BUILD_QUEUE_FIELDNAMES,
        MINIMUM_INDEPENDENT_SEEDS,
        _safe_identifier,
        _split_joined_values,
    )
except ModuleNotFoundError:  # pragma: no cover - fallback when run as scripts/*.py
    from p05_plan_family_profiles import (
        BUILD_QUEUE_FIELDNAMES,
        MINIMUM_INDEPENDENT_SEEDS,
        _safe_identifier,
        _split_joined_values,
    )


OUTPUT_COMMAND_MANIFEST_FILENAME = "p05_family_profile_command_manifest.tsv"
OUTPUT_COMMAND_SUMMARY_FILENAME = "p05_family_profile_command_summary.tsv"
COMMAND_FIELDNAMES = (
    "family_category",
    "seed_bundle_path",
    "alignment_path",
    "hmm_path",
    "alignment_command",
    "hmmbuild_command",
    "command_status",
    "notes",
)
SUMMARY_FIELDNAMES = ("kind", "name", "count")
COMMAND_STATUS = "planned_not_run"


def load_scaffold_queue(path: Path) -> list[dict[str, str]]:
    """Load and minimally validate the saved P05 scaffold queue."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing_columns = [field for field in BUILD_QUEUE_FIELDNAMES if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

        rows: list[dict[str, str]] = []
        seen_families: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            for field in (
                "family_category",
                "seed_bundle_path",
                "seed_row_count",
                "qualifying_seed_row_count",
                "independent_qualifying_accession_count",
                "qualifying_source_accessions",
                "bundled_sequence_count",
            ):
                if not normalized.get(field):
                    raise ValueError(f"{path}:{line_number} missing required value for {field}")
            family_category = normalized["family_category"]
            if family_category in seen_families:
                raise ValueError(f"{path}:{line_number} has duplicate family_category {family_category!r}")
            seen_families.add(family_category)
            seed_row_count = _parse_positive_int(normalized["seed_row_count"], path=path, line_number=line_number, field="seed_row_count")
            qualifying_seed_row_count = _parse_positive_int(
                normalized["qualifying_seed_row_count"], path=path, line_number=line_number, field="qualifying_seed_row_count"
            )
            independent_qualifying_accession_count = _parse_positive_int(
                normalized["independent_qualifying_accession_count"],
                path=path,
                line_number=line_number,
                field="independent_qualifying_accession_count",
            )
            bundled_sequence_count = _parse_positive_int(
                normalized["bundled_sequence_count"], path=path, line_number=line_number, field="bundled_sequence_count"
            )
            qualifying_source_accessions = _split_joined_values(normalized["qualifying_source_accessions"])
            if independent_qualifying_accession_count < MINIMUM_INDEPENDENT_SEEDS:
                raise ValueError(
                    f"{path}:{line_number} has only {independent_qualifying_accession_count} independent qualifying accessions; "
                    f"at least {MINIMUM_INDEPENDENT_SEEDS} are required"
                )
            if len(qualifying_source_accessions) != independent_qualifying_accession_count:
                raise ValueError(
                    f"{path}:{line_number} qualifying_source_accessions does not match independent_qualifying_accession_count"
                )
            if qualifying_seed_row_count < independent_qualifying_accession_count:
                raise ValueError(
                    f"{path}:{line_number} qualifying_seed_row_count cannot be smaller than independent_qualifying_accession_count"
                )
            if seed_row_count < qualifying_seed_row_count:
                raise ValueError(f"{path}:{line_number} seed_row_count cannot be smaller than qualifying_seed_row_count")
            if bundled_sequence_count != qualifying_seed_row_count:
                raise ValueError(
                    f"{path}:{line_number} bundled_sequence_count must match qualifying_seed_row_count"
                )
            rows.append(normalized)
    return rows


def build_family_profile_commands(queue_path: Path, outdir: Path) -> dict[str, Path]:
    """Write deterministic command and summary manifests from a scaffold TSV."""

    queue_rows = load_scaffold_queue(queue_path)
    alignment_dir = outdir.parent / "alignments"
    hmm_dir = outdir.parent / "hmms"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    hmm_dir.mkdir(parents=True, exist_ok=True)
    command_rows: list[dict[str, str]] = []
    for queue_row in sorted(queue_rows, key=lambda row: row["family_category"]):
        family_category = queue_row["family_category"]
        identifier = _safe_identifier(family_category)
        seed_bundle_path = _posix_path(queue_row["seed_bundle_path"])
        alignment_path = (alignment_dir / f"{identifier}.aligned.faa").as_posix()
        hmm_path = (hmm_dir / f"{identifier}.hmm").as_posix()
        command_rows.append(
            {
                "family_category": family_category,
                "seed_bundle_path": seed_bundle_path,
                "alignment_path": alignment_path,
                "hmm_path": hmm_path,
                "alignment_command": (
                    "mafft --localpair --maxiterate 1000 --inputorder "
                    f"{_shell_quote(seed_bundle_path)} > {_shell_quote(alignment_path)}"
                ),
                "hmmbuild_command": (
                    f"hmmbuild --amino {_shell_quote(hmm_path)} "
                    f"{_shell_quote(alignment_path)}"
                ),
                "command_status": COMMAND_STATUS,
                "notes": (
                    "Prepared from the saved P05 scaffold queue; MAFFT and hmmbuild were not run. "
                    "Calibrate against close non-target hydrolases before downstream use."
                ),
            }
        )

    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / OUTPUT_COMMAND_MANIFEST_FILENAME
    summary_path = outdir / OUTPUT_COMMAND_SUMMARY_FILENAME
    write_tsv(manifest_path, command_rows, COMMAND_FIELDNAMES)
    write_tsv(summary_path, summarize_command_manifest(queue_rows, command_rows), SUMMARY_FIELDNAMES)
    return {"manifest": manifest_path, "summary": summary_path}


def summarize_command_manifest(
    queue_rows: list[dict[str, str]], command_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    """Return a compact summary of command-manifest preparation."""

    return [
        {"kind": "total", "name": "scaffold_queue_rows", "count": str(len(queue_rows))},
        {"kind": "total", "name": "eligible_families", "count": str(len(command_rows))},
        {"kind": "total", "name": "command_manifest_rows", "count": str(len(command_rows))},
        {
            "kind": "command_status",
            "name": COMMAND_STATUS,
            "count": str(sum(row["command_status"] == COMMAND_STATUS for row in command_rows)),
        },
    ]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _shell_quote(value: str) -> str:
    """Quote a path for the POSIX shell used by the T141 runner."""

    return "'" + value.replace("'", "'\"'\"'") + "'"


def _posix_path(value: str) -> str:
    return value.replace("\\", "/")


def _parse_positive_int(value: str, *, path: Path, line_number: int, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"{path}:{line_number} {field} must be an integer") from error
    if parsed < 0:
        raise ValueError(f"{path}:{line_number} {field} must not be negative")
    return parsed


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare P05 MAFFT and HMMER command manifests.")
    parser.add_argument(
        "--scaffold-queue",
        type=Path,
        default=Path("04_family_profiles/manifests/p05_family_hmm_build_scaffold_queue.tsv"),
        help="Saved P05 scaffold queue TSV",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("04_family_profiles/manifests"),
        help="Output directory for the ignored command manifest and summary",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    outputs = build_family_profile_commands(args.scaffold_queue, args.outdir)
    print(f"Command manifest written: {outputs['manifest']}")
    print(f"Command summary written: {outputs['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
