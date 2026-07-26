"""Plan and scaffold P05 family-profile work from the normalized P04 reference manifest.

This stage keeps the biological boundary explicit:
custom HMMs are only planned for families with at least three independent
qualifying accessions. Bacterial rows still use direct/literature-backed seed
evidence, while archaeal rows may also count annotation-supported E3 records
when they are explicitly admitted under the archaeal library policy. Smaller
families remain anchor sets and do not get overfit HMMs from undersampled
evidence.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

try:  # pragma: no cover - direct script execution from repo root
    from scripts.p04_build_reference_library import load_reference_manifest
except ModuleNotFoundError:  # pragma: no cover - fallback when run as scripts/*.py
    from p04_build_reference_library import load_reference_manifest


MINIMUM_INDEPENDENT_SEEDS = 3
QUALIFYING_EVIDENCE_LEVELS_BY_DOMAIN = {
    "Bacteria": {"E1", "E2"},
    "Archaea": {"E1", "E2", "E3"},
}
DEFAULT_FAMILY_CLASSIFICATION_PATH = Path("04_family_profiles/manifests/p05_family_keep_now.tsv")
OUTPUT_PLAN_FILENAME = "p05_family_profile_plan.tsv"
OUTPUT_HMM_QUEUE_FILENAME = "p05_family_hmm_build_queue.tsv"
OUTPUT_ANCHOR_QUEUE_FILENAME = "p05_family_anchor_set_queue.tsv"
OUTPUT_SUMMARY_FILENAME = "p05_family_profile_summary.tsv"
OUTPUT_BUILD_QUEUE_FILENAME = "p05_family_hmm_build_scaffold_queue.tsv"
OUTPUT_BUILD_SUMMARY_FILENAME = "p05_family_hmm_build_scaffold_summary.tsv"
DEFAULT_SEED_BUNDLE_DIRNAME = "seed_bundles"
FASTA_LINE_WIDTH = 80

PLAN_FIELDNAMES = (
    "family_category",
    "taxonomic_domain",
    "reference_library",
    "seed_row_count",
    "qualifying_seed_row_count",
    "independent_qualifying_accession_count",
    "eligible_for_hmm",
    "profile_strategy",
    "evidence_levels",
    "qualifying_seed_ids",
    "qualifying_source_accessions",
    "recommended_next_step",
    "calibration_control_panel",
    "notes",
)

BUILD_QUEUE_FIELDNAMES = (
    "family_category",
    "taxonomic_domain",
    "reference_library",
    "seed_row_count",
    "qualifying_seed_row_count",
    "independent_qualifying_accession_count",
    "qualifying_seed_ids",
    "qualifying_source_accessions",
    "seed_bundle_path",
    "bundled_sequence_count",
    "alignment_tool",
    "alignment_mode",
    "hmm_build_tool",
    "calibration_search_tool",
    "calibration_control_panel",
    "source_manifest_path",
    "source_plan_path",
    "notes",
)

SUMMARY_FIELDNAMES = ("kind", "name", "count")


def plan_family_profiles(
    rows: list[dict[str, str]],
    minimum_independent_seeds: int = MINIMUM_INDEPENDENT_SEEDS,
) -> list[dict[str, str]]:
    """Convert reference-library rows into per-family profile plans."""

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["family_category"]].append(row)

    plan_rows: list[dict[str, str]] = []
    for family_category in sorted(grouped):
        family_rows = grouped[family_category]
        domains = sorted({row["taxonomic_domain"] for row in family_rows})
        if len(domains) != 1:
            raise ValueError(
                f"family_category {family_category!r} spans multiple taxonomic domains: {', '.join(domains)}"
            )

        reference_libraries = sorted({row["reference_library"] for row in family_rows})
        qualifying_evidence_levels = QUALIFYING_EVIDENCE_LEVELS_BY_DOMAIN[domains[0]]
        qualifying_rows = [row for row in family_rows if row["evidence_level"] in qualifying_evidence_levels]
        qualifying_seed_ids = sorted({row["seed_id"] for row in qualifying_rows})
        qualifying_accessions = sorted({row["source_accession"] for row in qualifying_rows})
        evidence_levels = _unique_join(row["evidence_level"] for row in family_rows)
        strategy = "build_hmm" if len(qualifying_accessions) >= minimum_independent_seeds else "anchor_set"

        if strategy == "build_hmm":
            next_step = "align_with_mafft_then_build_hmm_and_calibrate_with_close_non_target_hydrolases"
            calibration_panel = "close_non_target_hydrolases"
            notes = (
                f"{len(qualifying_accessions)} independent qualifying accessions meet the minimum "
                f"of {minimum_independent_seeds}; build a custom family HMM only after calibration."
            )
        else:
            next_step = "keep_as_anchor_set_and_expand_independent_experimental_seeds"
            calibration_panel = "not_applicable_yet"
            notes = (
                f"{len(qualifying_accessions)} independent qualifying accessions are below the "
                f"minimum of {minimum_independent_seeds}; do not build an overfit HMM."
            )

        plan_rows.append(
            {
                "family_category": family_category,
                "taxonomic_domain": domains[0],
                "reference_library": _unique_join(reference_libraries),
                "seed_row_count": str(len(family_rows)),
                "qualifying_seed_row_count": str(len(qualifying_rows)),
                "independent_qualifying_accession_count": str(len(qualifying_accessions)),
                "eligible_for_hmm": "yes" if strategy == "build_hmm" else "no",
                "profile_strategy": strategy,
                "evidence_levels": evidence_levels,
                "qualifying_seed_ids": _unique_join(qualifying_seed_ids),
                "qualifying_source_accessions": _unique_join(qualifying_accessions),
                "recommended_next_step": next_step,
                "calibration_control_panel": calibration_panel,
                "notes": notes,
            }
        )

    return plan_rows


def load_active_family_categories(path: Path) -> set[str]:
    """Load the family categories marked keep_now in a classification TSV."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing_columns = [field for field in ("family_category", "priority_status") if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

        active_categories: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            family_category = normalized.get("family_category", "")
            priority_status = normalized.get("priority_status", "")
            if not family_category:
                raise ValueError(f"{path}:{line_number} missing required value for family_category")
            if not priority_status:
                raise ValueError(f"{path}:{line_number} missing required value for priority_status")
            if priority_status not in {"keep_now", "deferred"}:
                raise ValueError(
                    f"{path}:{line_number} priority_status must be keep_now or deferred"
                )
            if priority_status == "keep_now":
                active_categories.add(family_category)

    if not active_categories:
        raise ValueError(f"{path} does not contain any keep_now family categories")
    return active_categories


def summarize_family_profile_plan(plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize the family-profile plan in a compact TSV."""

    summary: list[dict[str, str]] = []
    strategy_counts = Counter(row["profile_strategy"] for row in plan_rows)
    domain_counts = Counter(row["taxonomic_domain"] for row in plan_rows)
    summary.append({"kind": "total", "name": "families", "count": str(len(plan_rows))})
    summary.append(
        {
            "kind": "total",
            "name": "qualifying_seed_rows",
            "count": str(sum(int(row["qualifying_seed_row_count"]) for row in plan_rows)),
        }
    )
    summary.append(
        {
            "kind": "total",
            "name": "independent_qualifying_accessions",
            "count": str(sum(int(row["independent_qualifying_accession_count"]) for row in plan_rows)),
        }
    )
    for strategy in ("build_hmm", "anchor_set"):
        summary.append(
            {
                "kind": "profile_strategy",
                "name": strategy,
                "count": str(strategy_counts.get(strategy, 0)),
            }
        )
    for domain in ("Bacteria", "Archaea"):
        summary.append(
            {
                "kind": "taxonomic_domain",
                "name": domain,
                "count": str(domain_counts.get(domain, 0)),
            }
        )
    return summary


def load_family_profile_plan(path: Path) -> list[dict[str, str]]:
    """Load the planned family-profile TSV written by the planner."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing_columns = [field for field in PLAN_FIELDNAMES if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            validate_plan_row(normalized, path=path, line_number=line_number)
            rows.append(normalized)
    return rows


def validate_plan_row(row: dict[str, str], path: Path, line_number: int) -> None:
    """Validate a single P05 plan row."""

    for field in PLAN_FIELDNAMES:
        if not row.get(field):
            raise ValueError(f"{path}:{line_number} missing required value for {field}")

    if row["eligible_for_hmm"] not in {"yes", "no"}:
        raise ValueError(f"{path}:{line_number} eligible_for_hmm must be yes or no")

    if row["profile_strategy"] not in {"build_hmm", "anchor_set"}:
        raise ValueError(
            f"{path}:{line_number} profile_strategy must be build_hmm or anchor_set"
        )

    if row["profile_strategy"] == "build_hmm" and row["eligible_for_hmm"] != "yes":
        raise ValueError(
            f"{path}:{line_number} build_hmm rows must have eligible_for_hmm set to yes"
        )

    if row["profile_strategy"] == "anchor_set" and row["eligible_for_hmm"] != "no":
        raise ValueError(
            f"{path}:{line_number} anchor_set rows must have eligible_for_hmm set to no"
        )

    try:
        int(row["seed_row_count"])
        int(row["qualifying_seed_row_count"])
        int(row["independent_qualifying_accession_count"])
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} count columns must contain integers"
        ) from error


def build_family_profile_queue(
    manifest_rows: list[dict[str, str]],
    plan_rows: list[dict[str, str]],
    minimum_independent_seeds: int = MINIMUM_INDEPENDENT_SEEDS,
    manifest_path: Path | None = None,
    plan_path: Path | None = None,
) -> list[dict[str, str]]:
    """Build the MAFFT/HMMER queue from manifest rows and a saved P05 plan."""

    manifest_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        manifest_by_family[row["family_category"]].append(row)

    queue_rows: list[dict[str, str]] = []
    for plan_row in sorted(plan_rows, key=lambda row: row["family_category"]):
        if plan_row["profile_strategy"] != "build_hmm":
            continue

        family_category = plan_row["family_category"]
        family_rows = manifest_by_family.get(family_category, [])
        qualifying_evidence_levels = QUALIFYING_EVIDENCE_LEVELS_BY_DOMAIN[plan_row["taxonomic_domain"]]
        qualifying_rows = [row for row in family_rows if row["evidence_level"] in qualifying_evidence_levels]
        manifest_accessions = sorted({row["source_accession"] for row in qualifying_rows})
        plan_accessions = sorted(_split_joined_values(plan_row["qualifying_source_accessions"]))

        if manifest_accessions != plan_accessions:
            raise ValueError(
                f"family_category {family_category!r} plan and manifest disagree on qualifying source accessions"
            )

        if len(manifest_accessions) < minimum_independent_seeds:
            raise ValueError(
                f"family_category {family_category!r} has only {len(manifest_accessions)} independent qualifying accessions; "
                f"at least {minimum_independent_seeds} are required"
            )

        queue_rows.append(
            {
                "family_category": family_category,
                "taxonomic_domain": plan_row["taxonomic_domain"],
                "reference_library": plan_row["reference_library"],
                "seed_row_count": str(len(family_rows)),
                "qualifying_seed_row_count": str(len(qualifying_rows)),
                "independent_qualifying_accession_count": str(len(manifest_accessions)),
                "qualifying_seed_ids": _unique_join(row["seed_id"] for row in qualifying_rows),
                "qualifying_source_accessions": _unique_join(manifest_accessions),
                "seed_bundle_path": "",
                "bundled_sequence_count": "0",
                "alignment_tool": "MAFFT",
                "alignment_mode": "L-INS-i",
                "hmm_build_tool": "hmmbuild",
                "calibration_search_tool": "hmmsearch",
                "calibration_control_panel": "close_non_target_hydrolases",
                "source_manifest_path": _as_posix_path(manifest_path),
                "source_plan_path": _as_posix_path(plan_path),
                "notes": (
                    f"Three-independent-accession threshold met; align with MAFFT L-INS-i and build "
                    f"the family HMM with hmmbuild before calibration."
                ),
            }
        )

    return queue_rows


def summarize_family_profile_queue(
    manifest_rows: list[dict[str, str]],
    plan_rows: list[dict[str, str]],
    queue_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Summarize the build queue in a compact TSV."""

    manifest_families = {row["family_category"] for row in manifest_rows}
    plan_families = {row["family_category"] for row in plan_rows}
    queue_families = {row["family_category"] for row in queue_rows}
    summary: list[dict[str, str]] = [
        {"kind": "total", "name": "families_in_manifest", "count": str(len(manifest_families))},
        {"kind": "total", "name": "families_in_plan", "count": str(len(plan_families))},
        {"kind": "total", "name": "eligible_families", "count": str(len(queue_families))},
        {"kind": "total", "name": "ineligible_families", "count": str(len(plan_families - queue_families))},
        {"kind": "total", "name": "build_queue_rows", "count": str(len(queue_rows))},
        {
            "kind": "total",
            "name": "qualifying_seed_rows",
            "count": str(sum(int(row["qualifying_seed_row_count"]) for row in queue_rows)),
        },
        {
            "kind": "total",
            "name": "independent_qualifying_accessions",
            "count": str(sum(int(row["independent_qualifying_accession_count"]) for row in queue_rows)),
        },
    ]
    return summary


def build_family_profile_scaffold(
    manifest_path: Path,
    plan_path: Path,
    outdir: Path,
    minimum_independent_seeds: int = MINIMUM_INDEPENDENT_SEEDS,
    classification_path: Path | None = DEFAULT_FAMILY_CLASSIFICATION_PATH,
    bundle_dir: Path | None = None,
) -> dict[str, Path]:
    """Write the build queue and deterministic seed bundles for eligible families."""

    manifest_rows = load_reference_manifest(manifest_path)
    if not plan_path.exists():
        plan_outputs = build_family_profile_plan(
            manifest_path,
            outdir,
            minimum_independent_seeds=minimum_independent_seeds,
            classification_path=classification_path,
        )
        plan_path = plan_outputs["plan"]
    plan_rows = load_family_profile_plan(plan_path)
    queue_rows = build_family_profile_queue(
        manifest_rows,
        plan_rows,
        minimum_independent_seeds=minimum_independent_seeds,
        manifest_path=manifest_path,
        plan_path=plan_path,
    )
    bundle_dir = bundle_dir or outdir.parent / DEFAULT_SEED_BUNDLE_DIRNAME
    queue_rows = materialize_family_seed_bundles(
        manifest_rows,
        queue_rows,
        manifest_path=manifest_path,
        bundle_dir=bundle_dir,
    )
    summary_rows = summarize_family_profile_queue(manifest_rows, plan_rows, queue_rows)

    outdir.mkdir(parents=True, exist_ok=True)
    queue_path = outdir / OUTPUT_BUILD_QUEUE_FILENAME
    summary_path = outdir / OUTPUT_BUILD_SUMMARY_FILENAME

    write_tsv(queue_path, queue_rows, BUILD_QUEUE_FIELDNAMES)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)

    return {
        "queue": queue_path,
        "summary": summary_path,
    }


def materialize_family_seed_bundles(
    manifest_rows: list[dict[str, str]],
    queue_rows: list[dict[str, str]],
    manifest_path: Path,
    bundle_dir: Path,
) -> list[dict[str, str]]:
    """Write one deterministic unaligned FASTA bundle per eligible family."""

    if not queue_rows:
        return queue_rows

    manifest_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        qualifying_evidence_levels = QUALIFYING_EVIDENCE_LEVELS_BY_DOMAIN[row["taxonomic_domain"]]
        if row["evidence_level"] in qualifying_evidence_levels:
            manifest_by_family[row["family_category"]].append(row)

    bundle_dir.mkdir(parents=True, exist_ok=True)
    materialized_rows: list[dict[str, str]] = []
    for queue_row in queue_rows:
        family_category = queue_row["family_category"]
        family_rows = sorted(
            manifest_by_family.get(family_category, []),
            key=lambda row: (row["seed_id"], row["source_accession"], row["sequence_path"]),
        )
        records: list[tuple[str, str]] = []
        seen_headers: set[str] = set()
        for seed_row in family_rows:
            header = stable_seed_fasta_header(seed_row["seed_id"], seed_row["source_accession"])
            if header in seen_headers:
                raise ValueError(
                    f"family_category {family_category!r} has duplicate deterministic FASTA header {header!r}"
            )
            sequence_path = _resolve_manifest_sequence_path(manifest_path, seed_row["sequence_path"])
            sequence = _read_single_fasta_sequence(sequence_path)
            records.append((header, sequence))
            seen_headers.add(header)

        bundle_path = bundle_dir / f"{_safe_identifier(family_category)}.faa"
        _write_fasta(bundle_path, records)
        materialized_row = dict(queue_row)
        materialized_row["seed_bundle_path"] = _as_posix_path(bundle_path)
        materialized_row["bundled_sequence_count"] = str(len(records))
        materialized_rows.append(materialized_row)

    return materialized_rows


def stable_seed_fasta_header(seed_id: str, source_accession: str) -> str:
    """Return the stable FASTA identifier used by P05 seed bundles."""

    return f"{_safe_identifier(seed_id)}|{_safe_identifier(source_accession)}"


def _resolve_manifest_sequence_path(manifest_path: Path, sequence_path: str) -> Path:
    path = Path(sequence_path)
    if path.is_absolute():
        return path

    candidate_paths = [manifest_path.parent / path]
    repo_root = manifest_path.parent.parent
    if repo_root != manifest_path.parent:
        candidate_paths.append(repo_root / path)

    for candidate in candidate_paths:
        if candidate.exists():
            return candidate

    return candidate_paths[0]


def _as_posix_path(path: Path | None) -> str:
    if path is None:
        return ""
    return path.as_posix()


def _read_single_fasta_sequence(path: Path) -> str:
    header_seen = False
    sequence_parts: list[str] = []
    record_count = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            record_count += 1
            if record_count > 1:
                raise ValueError(f"{path} must contain exactly one FASTA record")
            header_seen = bool(line[1:].strip())
            if not header_seen:
                raise ValueError(f"{path} has an empty FASTA header")
            continue
        if not header_seen:
            raise ValueError(f"{path} contains sequence data before a FASTA header")
        sequence_parts.append("".join(line.split()))

    if record_count != 1:
        raise ValueError(f"{path} must contain exactly one FASTA record")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"{path} has an empty FASTA sequence")
    return sequence


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    chunks: list[str] = []
    for header, sequence in records:
        chunks.append(f">{header}\n")
        chunks.extend(
            f"{sequence[start:start + FASTA_LINE_WIDTH]}\n"
            for start in range(0, len(sequence), FASTA_LINE_WIDTH)
        )
    path.write_text("".join(chunks), encoding="utf-8", newline="\n")


def _safe_identifier(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._:-]+", "_", value.strip())
    return sanitized or "unnamed"


def build_family_profile_plan(
    manifest_path: Path,
    outdir: Path,
    minimum_independent_seeds: int = MINIMUM_INDEPENDENT_SEEDS,
    classification_path: Path | None = DEFAULT_FAMILY_CLASSIFICATION_PATH,
) -> dict[str, Path]:
    """Read the normalized reference manifest and write P05 planning TSVs."""

    rows = load_reference_manifest(manifest_path)
    if classification_path is not None:
        active_family_categories = load_active_family_categories(classification_path)
        manifest_family_categories = {row["family_category"] for row in rows}
        missing_active_families = sorted(active_family_categories - manifest_family_categories)
        if missing_active_families:
            raise ValueError(
                "keep_now family categories not present in the reference manifest: "
                + ", ".join(missing_active_families)
            )
        rows = [row for row in rows if row["family_category"] in active_family_categories]
    plan_rows = plan_family_profiles(rows, minimum_independent_seeds=minimum_independent_seeds)
    summary_rows = summarize_family_profile_plan(plan_rows)
    hmm_ready_rows = [row for row in plan_rows if row["eligible_for_hmm"] == "yes"]
    anchor_rows = [row for row in plan_rows if row["eligible_for_hmm"] == "no"]

    outdir.mkdir(parents=True, exist_ok=True)
    plan_path = outdir / OUTPUT_PLAN_FILENAME
    hmm_ready_path = outdir / OUTPUT_HMM_QUEUE_FILENAME
    anchor_path = outdir / OUTPUT_ANCHOR_QUEUE_FILENAME
    summary_path = outdir / OUTPUT_SUMMARY_FILENAME

    write_tsv(plan_path, plan_rows, PLAN_FIELDNAMES)
    write_tsv(hmm_ready_path, hmm_ready_rows, PLAN_FIELDNAMES)
    write_tsv(anchor_path, anchor_rows, PLAN_FIELDNAMES)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)

    return {
        "plan": plan_path,
        "hmm_ready": hmm_ready_path,
        "anchor_sets": anchor_path,
        "summary": summary_path,
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plan P05 family-profile work from the P04 reference manifest.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("01_reference_library/reference_library.normalized.tsv"),
        help="Normalized P04 reference manifest TSV",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("04_family_profiles/manifests"),
        help="Output directory for P05 TSVs",
    )
    parser.add_argument(
        "--plan-path",
        type=Path,
        default=Path("04_family_profiles/manifests/p05_family_profile_plan.tsv"),
        help="Existing P05 plan TSV to use for the build scaffold",
    )
    parser.add_argument(
        "--minimum-independent-seeds",
        type=int,
        default=MINIMUM_INDEPENDENT_SEEDS,
        help="Minimum distinct qualifying source accessions required to plan an HMM",
    )
    parser.add_argument(
        "--family-classification",
        type=Path,
        default=DEFAULT_FAMILY_CLASSIFICATION_PATH,
        help="TSV listing the family categories to keep active for P05 planning",
    )
    parser.add_argument(
        "--build-scaffold",
        action="store_true",
        help="Build the MAFFT/HMMER scaffold queue from the saved P05 plan",
    )
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        default=None,
        help="Directory for per-family seed FASTA bundles (defaults beside the scaffold manifest directory)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.build_scaffold:
        outputs = build_family_profile_scaffold(
            args.manifest,
            args.plan_path,
            args.outdir,
            minimum_independent_seeds=args.minimum_independent_seeds,
            classification_path=args.family_classification,
            bundle_dir=args.bundle_dir,
        )
        print(f"Build queue written: {outputs['queue']}")
        print(f"Build summary written: {outputs['summary']}")
    else:
        outputs = build_family_profile_plan(
            args.manifest,
            args.outdir,
            minimum_independent_seeds=args.minimum_independent_seeds,
            classification_path=args.family_classification,
        )
        print(f"Family-profile plan written: {outputs['plan']}")
        print(f"HMM-ready queue written: {outputs['hmm_ready']}")
        print(f"Anchor-set queue written: {outputs['anchor_sets']}")
        print(f"Summary written: {outputs['summary']}")
    return 0


def _unique_join(values: list[str] | tuple[str, ...] | set[str] | object) -> str:
    unique = sorted({value for value in values if value})
    return ";".join(unique)


def _split_joined_values(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


if __name__ == "__main__":
    raise SystemExit(main())
