"""Audit P06 candidate-table reasonableness before P07 interpretation.

The audit summarizes P06 HMMER candidate rows into compact, reproducible
counts. It checks scale and overlap, but deliberately does not convert HMM
hits into PHB/PHA degradation phenotypes.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


OUTPUT_SUMMARY_FILENAME = "p06_candidate_reasonableness_summary.tsv"
OUTPUT_FAMILY_TIER_FILENAME = "p06_family_tier_reasonableness.tsv"
OUTPUT_OVERLAP_FILENAME = "p06_high_confidence_overlap_targets.tsv"
OUTPUT_REPORT_FILENAME = "P06_REASONABLENESS_AUDIT.md"
REQUIRED_CANDIDATE_FIELDS = (
    "family_category",
    "proteome_shard",
    "target_id",
    "full_sequence_score",
    "tier",
)
SUMMARY_FIELDNAMES = ("kind", "name", "value")
FAMILY_TIER_FIELDNAMES = (
    "family_category",
    "tier",
    "candidate_rows",
    "unique_targets",
    "proteome_shards",
    "max_full_sequence_score",
)
OVERLAP_FIELDNAMES = (
    "target_id",
    "high_confidence_family_count",
    "family_categories",
    "proteome_shards",
    "p06_candidate_rows",
    "notes",
)


def audit_p06_candidates(
    candidate_table: Path,
    outdir: Path,
    *,
    total_predicted_genes: int | None = None,
    total_genomes: int | None = None,
) -> dict[str, Path]:
    """Write compact P06 reasonableness summaries."""

    rows = load_candidate_rows(candidate_table)
    outdir.mkdir(parents=True, exist_ok=True)

    family_tier_rows = build_family_tier_rows(rows)
    overlap_rows = build_high_confidence_overlap_rows(rows)
    summary_rows = build_summary_rows(
        rows,
        overlap_rows,
        total_predicted_genes=total_predicted_genes,
        total_genomes=total_genomes,
    )
    report_text = build_markdown_report(
        rows,
        summary_rows,
        family_tier_rows,
        overlap_rows,
        candidate_table=candidate_table,
        total_predicted_genes=total_predicted_genes,
        total_genomes=total_genomes,
    )

    summary_path = outdir / OUTPUT_SUMMARY_FILENAME
    family_tier_path = outdir / OUTPUT_FAMILY_TIER_FILENAME
    overlap_path = outdir / OUTPUT_OVERLAP_FILENAME
    report_path = outdir / OUTPUT_REPORT_FILENAME
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)
    write_tsv(family_tier_path, family_tier_rows, FAMILY_TIER_FIELDNAMES)
    write_tsv(overlap_path, overlap_rows, OVERLAP_FIELDNAMES)
    report_path.write_text(report_text, encoding="utf-8")

    return {
        "summary": summary_path,
        "family_tier_summary": family_tier_path,
        "overlap_summary": overlap_path,
        "markdown_report": report_path,
    }


def load_candidate_rows(candidate_table: Path) -> list[dict[str, str]]:
    """Load P06 candidate rows and validate the minimum audit schema."""

    with candidate_table.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{candidate_table} is missing a tabular header")
        missing = [field for field in REQUIRED_CANDIDATE_FIELDS if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"{candidate_table} is missing required columns: {', '.join(missing)}")

        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            for field_name in REQUIRED_CANDIDATE_FIELDS:
                if not normalized.get(field_name):
                    raise ValueError(f"{candidate_table}:{line_number} missing required value for {field_name}")
            score = _parse_score(normalized["full_sequence_score"], candidate_table, line_number)
            normalized["full_sequence_score"] = f"{score:.1f}"
            rows.append(normalized)
    if not rows:
        raise ValueError(f"{candidate_table} has no candidate rows to audit")
    return rows


def build_family_tier_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize row, target, and shard counts by family and tier."""

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["family_category"], row["tier"])].append(row)

    output_rows: list[dict[str, str]] = []
    for family, tier in sorted(grouped):
        group_rows = grouped[(family, tier)]
        output_rows.append(
            {
                "family_category": family,
                "tier": tier,
                "candidate_rows": str(len(group_rows)),
                "unique_targets": str(len({row["target_id"] for row in group_rows})),
                "proteome_shards": str(len({row["proteome_shard"] for row in group_rows})),
                "max_full_sequence_score": f"{max(float(row['full_sequence_score']) for row in group_rows):.1f}",
            }
        )
    return output_rows


def build_high_confidence_overlap_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Find High-confidence targets assigned to more than one P06 model."""

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["tier"] == "High-confidence":
            grouped[row["target_id"]].append(row)

    overlap_rows: list[dict[str, str]] = []
    for target_id in sorted(grouped):
        group_rows = grouped[target_id]
        families = sorted({row["family_category"] for row in group_rows})
        if len(families) <= 1:
            continue
        overlap_rows.append(
            {
                "target_id": target_id,
                "high_confidence_family_count": str(len(families)),
                "family_categories": _unique_join(families),
                "proteome_shards": _unique_join(sorted({row["proteome_shard"] for row in group_rows})),
                "p06_candidate_rows": str(len(group_rows)),
                "notes": (
                    "High-confidence overlap is a sequence-classification review flag; "
                    "P07/P08 must resolve architecture and phylogeny before interpretation."
                ),
            }
        )
    return overlap_rows


def build_summary_rows(
    rows: list[dict[str, str]],
    overlap_rows: list[dict[str, str]],
    *,
    total_predicted_genes: int | None,
    total_genomes: int | None,
) -> list[dict[str, str]]:
    """Build compact overall counts and rates."""

    tiers = sorted({row["tier"] for row in rows})
    summary_rows = [
        {"kind": "total", "name": "candidate_rows", "value": str(len(rows))},
        {"kind": "total", "name": "unique_candidate_targets", "value": str(len({row["target_id"] for row in rows}))},
        {"kind": "total", "name": "proteome_shards", "value": str(len({row["proteome_shard"] for row in rows}))},
        {"kind": "total", "name": "families", "value": str(len({row["family_category"] for row in rows}))},
        {"kind": "total", "name": "high_confidence_overlap_targets", "value": str(len(overlap_rows))},
    ]
    for tier in tiers:
        tier_rows = [row for row in rows if row["tier"] == tier]
        summary_rows.append({"kind": "tier", "name": f"{tier}_rows", "value": str(len(tier_rows))})
        summary_rows.append(
            {
                "kind": "tier",
                "name": f"{tier}_unique_targets",
                "value": str(len({row["target_id"] for row in tier_rows})),
            }
        )
    high_rows = [row for row in rows if row["tier"] == "High-confidence"]
    high_unique_targets = len({row["target_id"] for row in high_rows})
    if total_predicted_genes is not None:
        _validate_positive_integer(total_predicted_genes, "total_predicted_genes")
        summary_rows.append(
            {
                "kind": "rate",
                "name": "High-confidence_rows_per_predicted_gene",
                "value": _format_rate(len(high_rows), total_predicted_genes),
            }
        )
        summary_rows.append(
            {
                "kind": "rate",
                "name": "High-confidence_unique_targets_per_predicted_gene",
                "value": _format_rate(high_unique_targets, total_predicted_genes),
            }
        )
    if total_genomes is not None:
        _validate_positive_integer(total_genomes, "total_genomes")
        summary_rows.append(
            {
                "kind": "rate",
                "name": "High-confidence_rows_per_genome",
                "value": _format_rate(len(high_rows), total_genomes),
            }
        )
        summary_rows.append(
            {
                "kind": "rate",
                "name": "High-confidence_unique_targets_per_genome",
                "value": _format_rate(high_unique_targets, total_genomes),
            }
        )
    return summary_rows


def build_markdown_report(
    rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    family_tier_rows: list[dict[str, str]],
    overlap_rows: list[dict[str, str]],
    *,
    candidate_table: Path,
    total_predicted_genes: int | None,
    total_genomes: int | None,
) -> str:
    """Write a short human-readable audit report."""

    generated_at = datetime.now(UTC).replace(microsecond=0).isoformat()
    summary_lookup = {(row["kind"], row["name"]): row["value"] for row in summary_rows}
    lines = [
        "# P06 candidate reasonableness audit",
        "",
        f"**Generated at UTC:** {generated_at}",
        f"**Candidate table:** `{candidate_table.as_posix()}`",
        "",
        "## Scope",
        "",
        (
            "This audit checks P06 HMMER candidate-table scale, de-duplicated target counts, "
            "and cross-model overlap. HMM hits are sequence evidence and not phenotype proof."
        ),
        "",
        "## Overall counts",
        "",
        f"- Candidate rows: `{summary_lookup[('total', 'candidate_rows')]}`",
        f"- Unique target ids: `{summary_lookup[('total', 'unique_candidate_targets')]}`",
        f"- Families: `{summary_lookup[('total', 'families')]}`",
        f"- Proteome shards represented: `{summary_lookup[('total', 'proteome_shards')]}`",
        f"- High-confidence overlap targets: `{summary_lookup[('total', 'high_confidence_overlap_targets')]}`",
    ]
    if total_predicted_genes is not None:
        lines.append(f"- P03 predicted genes denominator: `{total_predicted_genes}`")
        lines.append(
            "- High-confidence rows / predicted gene: "
            f"`{summary_lookup[('rate', 'High-confidence_rows_per_predicted_gene')]}`"
        )
    if total_genomes is not None:
        lines.append(f"- GTDB representative genome denominator: `{total_genomes}`")
        lines.append(
            "- High-confidence unique targets / genome: "
            f"`{summary_lookup[('rate', 'High-confidence_unique_targets_per_genome')]}`"
        )
    lines.extend(
        [
            "",
            "## Family and tier counts",
            "",
            "| Family | Tier | Rows | Unique targets | Shards | Max score |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in family_tier_rows:
        lines.append(
            "| "
            f"{row['family_category']} | {row['tier']} | {row['candidate_rows']} | "
            f"{row['unique_targets']} | {row['proteome_shards']} | {row['max_full_sequence_score']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            (
                "Large High-confidence families should be carried into P07 domain architecture, localization, "
                "and later P08 taxonomy/phylogeny checks. Cross-model overlaps are review flags, not automatic "
                "false positives or phenotype assignments."
            ),
        ]
    )
    if overlap_rows:
        lines.extend(["", "## High-confidence overlap preview", ""])
        for row in overlap_rows[:20]:
            lines.append(
                f"- `{row['target_id']}`: {row['family_categories']} "
                f"across {row['proteome_shards']}"
            )
    return "\n".join(lines) + "\n"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _parse_score(value: str, candidate_table: Path, line_number: int) -> float:
    try:
        score = float(value)
    except ValueError as exc:
        raise ValueError(f"{candidate_table}:{line_number} has invalid full_sequence_score {value!r}") from exc
    if not math.isfinite(score):
        raise ValueError(f"{candidate_table}:{line_number} has non-finite full_sequence_score {value!r}")
    return score


def _validate_positive_integer(value: int, name: str) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")


def _format_rate(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.6f}"


def _unique_join(values: Iterable[str]) -> str:
    return ";".join(dict.fromkeys(values))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Audit P06 candidate-table reasonableness before P07.")
    parser.add_argument("--candidate-table", type=Path, default=Path("05_hmmer_scan/p06_hmmer_candidates.tsv"))
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("06_domain_annotation/p06_reasonableness"),
        help="Directory for compact P06 audit summaries",
    )
    parser.add_argument("--total-predicted-genes", type=int, default=None)
    parser.add_argument("--total-genomes", type=int, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    outputs = audit_p06_candidates(
        args.candidate_table,
        args.outdir,
        total_predicted_genes=args.total_predicted_genes,
        total_genomes=args.total_genomes,
    )
    print(f"Reasonableness summary written: {outputs['summary']}")
    print(f"Family/tier summary written: {outputs['family_tier_summary']}")
    print(f"Overlap summary written: {outputs['overlap_summary']}")
    print(f"Markdown report written: {outputs['markdown_report']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
