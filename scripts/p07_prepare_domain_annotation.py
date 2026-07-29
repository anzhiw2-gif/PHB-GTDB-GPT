"""Prepare P07 domain-architecture and localization annotation inputs.

P07 starts from P06 sequence-evidence candidates. This script extracts the
selected protein sequences from the GTDB R232 P03 proteomes, writes auditable
FASTA shards, and prepares command manifests for annotation tools. It does not
interpret HMM hits as PHB/PHA degradation phenotypes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable, Iterator, TextIO


OUTPUT_SEQUENCE_MANIFEST_FILENAME = "p07_candidate_sequence_manifest.tsv"
OUTPUT_COMMAND_MANIFEST_FILENAME = "p07_domain_annotation_command_manifest.tsv"
OUTPUT_SUMMARY_FILENAME = "p07_domain_annotation_summary.tsv"
OUTPUT_MISSING_FILENAME = "p07_missing_candidate_sequences.tsv"
INPUT_FASTA_DIRNAME = "input/fasta_shards"
INTERPRO_DIRNAME = "interpro"
LOCALIZATION_DIRNAME = "localization"
REVIEW_DIRNAME = "review"
COMMAND_STATUS = "planned_not_run"
DEFAULT_INCLUDE_TIERS = ("High-confidence",)
P06_CANDIDATE_REQUIRED_FIELDS = (
    "family_category",
    "proteome_shard",
    "target_id",
    "target_length",
    "full_sequence_score",
    "tier",
)
P06_SCAN_MANIFEST_REQUIRED_FIELDS = ("proteome_shard", "proteome_path")

SEQUENCE_MANIFEST_FIELDNAMES = (
    "p07_sequence_id",
    "proteome_shard",
    "target_id",
    "source_proteome_path",
    "target_length_from_p06",
    "sequence_length",
    "family_categories",
    "tiers",
    "max_full_sequence_score",
    "p06_candidate_rows",
    "fasta_shard",
    "candidate_table_path",
    "scan_manifest_path",
    "gtdb_release",
    "generated_at_utc",
    "notes",
)
COMMAND_MANIFEST_FIELDNAMES = (
    "tool",
    "fasta_shard",
    "sequence_count",
    "input_fasta",
    "output_path",
    "command",
    "command_status",
    "notes",
)
SUMMARY_FIELDNAMES = ("kind", "name", "count")
MISSING_FIELDNAMES = (
    "proteome_shard",
    "target_id",
    "expected_proteome_paths",
    "candidate_table_path",
    "scan_manifest_path",
    "notes",
)


@dataclass
class CandidateGroup:
    proteome_shard: str
    target_id: str
    families: set[str] = field(default_factory=set)
    tiers: set[str] = field(default_factory=set)
    target_lengths: set[str] = field(default_factory=set)
    scores: list[float] = field(default_factory=list)
    p06_candidate_rows: int = 0

    def add_row(self, row: dict[str, str]) -> None:
        self.families.add(row["family_category"])
        self.tiers.add(row["tier"])
        self.target_lengths.add(row["target_length"])
        self.p06_candidate_rows += 1
        try:
            score = float(row["full_sequence_score"])
        except ValueError as exc:
            raise ValueError(f"Invalid P06 full_sequence_score for {self.target_id!r}") from exc
        if not math.isfinite(score):
            raise ValueError(f"Non-finite P06 full_sequence_score for {self.target_id!r}")
        self.scores.append(score)

    @property
    def max_score(self) -> float:
        return max(self.scores) if self.scores else 0.0

    @property
    def p06_target_length(self) -> str:
        if len(self.target_lengths) != 1:
            values = ", ".join(sorted(self.target_lengths))
            raise ValueError(
                f"P06 target_length mismatch for {self.proteome_shard}/{self.target_id}: {values}"
            )
        return next(iter(self.target_lengths))


@dataclass
class FastaRecord:
    record_id: str
    header: str
    sequence: str
    source_path: Path


def prepare_p07_inputs(
    candidate_table: Path,
    scan_manifest: Path,
    outdir: Path,
    *,
    include_tiers: Iterable[str] = DEFAULT_INCLUDE_TIERS,
    sequences_per_shard: int = 500,
    gtdb_release: str = "GTDB Release 11 R232",
    interproscan_exe: str = "interproscan.sh",
    interproscan_cpu: int = 8,
    signalp_exe: str = "signalp6",
    signalp_mode: str = "fast",
    phobius_exe: str = "",
    allow_missing: bool = False,
) -> dict[str, Path]:
    """Extract selected P06 candidates and write P07 planning manifests."""

    if sequences_per_shard < 1:
        raise ValueError("sequences_per_shard must be at least 1")
    if interproscan_cpu < 1:
        raise ValueError("interproscan_cpu must be at least 1")

    selected_tiers = tuple(include_tiers)
    candidate_groups = load_candidate_groups(candidate_table, include_tiers=selected_tiers)
    shard_to_paths = load_scan_shard_paths(scan_manifest)
    _ensure_shards_have_paths(candidate_groups, shard_to_paths, scan_manifest)

    generated_at = datetime.now(UTC).replace(microsecond=0).isoformat()
    fasta_dir = outdir / INPUT_FASTA_DIRNAME
    manifest_dir = outdir / "manifests"
    review_dir = outdir / REVIEW_DIRNAME
    interpro_dir = outdir / INTERPRO_DIRNAME
    localization_dir = outdir / LOCALIZATION_DIRNAME
    for directory in (fasta_dir, manifest_dir, review_dir, interpro_dir, localization_dir):
        directory.mkdir(parents=True, exist_ok=True)

    found_records, missing_rows = extract_records(candidate_groups, shard_to_paths)
    if missing_rows and not allow_missing:
        missing_path = review_dir / OUTPUT_MISSING_FILENAME
        write_tsv(missing_path, missing_rows, MISSING_FIELDNAMES)
        raise ValueError(
            f"{len(missing_rows)} selected P06 candidate sequences were not found in the referenced proteomes; "
            f"details written to {missing_path}"
        )

    sequence_rows, shard_rows = write_fasta_shards(
        found_records,
        fasta_dir,
        candidate_table=candidate_table,
        scan_manifest=scan_manifest,
        gtdb_release=gtdb_release,
        generated_at_utc=generated_at,
        sequences_per_shard=sequences_per_shard,
    )
    command_rows = build_command_manifest_rows(
        shard_rows,
        outdir,
        interproscan_exe=interproscan_exe,
        interproscan_cpu=interproscan_cpu,
        signalp_exe=signalp_exe,
        signalp_mode=signalp_mode,
        phobius_exe=phobius_exe,
    )
    summary_rows = summarize_p07_plan(
        candidate_groups,
        sequence_rows,
        command_rows,
        selected_tiers=selected_tiers,
        missing_rows=missing_rows,
    )

    sequence_manifest_path = manifest_dir / OUTPUT_SEQUENCE_MANIFEST_FILENAME
    command_manifest_path = manifest_dir / OUTPUT_COMMAND_MANIFEST_FILENAME
    summary_path = manifest_dir / OUTPUT_SUMMARY_FILENAME
    write_tsv(sequence_manifest_path, sequence_rows, SEQUENCE_MANIFEST_FIELDNAMES)
    write_tsv(command_manifest_path, command_rows, COMMAND_MANIFEST_FIELDNAMES)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)
    if missing_rows:
        write_tsv(review_dir / OUTPUT_MISSING_FILENAME, missing_rows, MISSING_FIELDNAMES)

    return {
        "sequence_manifest": sequence_manifest_path,
        "command_manifest": command_manifest_path,
        "summary": summary_path,
    }


def load_candidate_groups(candidate_table: Path, *, include_tiers: Iterable[str]) -> dict[tuple[str, str], CandidateGroup]:
    """Load and de-duplicate selected P06 candidate rows by shard and target."""

    selected = {tier.strip() for tier in include_tiers if tier.strip()}
    if not selected:
        raise ValueError("At least one tier must be selected for P07 input preparation")

    with candidate_table.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_fields(candidate_table, reader.fieldnames, P06_CANDIDATE_REQUIRED_FIELDS)
        groups: dict[tuple[str, str], CandidateGroup] = {}
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            for field_name in P06_CANDIDATE_REQUIRED_FIELDS:
                if not normalized.get(field_name):
                    raise ValueError(f"{candidate_table}:{line_number} missing required value for {field_name}")
            if normalized["tier"] not in selected:
                continue
            key = (normalized["proteome_shard"], normalized["target_id"])
            group = groups.setdefault(
                key,
                CandidateGroup(proteome_shard=normalized["proteome_shard"], target_id=normalized["target_id"]),
            )
            group.add_row(normalized)
    if not groups:
        raise ValueError(f"{candidate_table} has no rows matching selected tiers: {', '.join(sorted(selected))}")
    return groups


def load_scan_shard_paths(scan_manifest: Path) -> dict[str, list[Path]]:
    """Map each P06 proteome shard or chunk to the exact proteome paths it scanned."""

    with scan_manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_fields(scan_manifest, reader.fieldnames, P06_SCAN_MANIFEST_REQUIRED_FIELDS)
        paths_by_shard: dict[str, tuple[str, ...]] = {}
        for line_number, row in enumerate(reader, start=2):
            shard = (row.get("proteome_shard") or "").strip()
            raw_paths = (row.get("proteome_path") or "").strip()
            if not shard:
                raise ValueError(f"{scan_manifest}:{line_number} missing required value for proteome_shard")
            if not raw_paths:
                raise ValueError(f"{scan_manifest}:{line_number} missing required value for proteome_path")
            paths = tuple(item.strip() for item in raw_paths.split(";") if item.strip())
            if shard in paths_by_shard and paths_by_shard[shard] != paths:
                raise ValueError(f"{scan_manifest}:{line_number} has inconsistent proteome_path for shard {shard!r}")
            paths_by_shard[shard] = paths
    return {shard: [Path(value) for value in paths] for shard, paths in paths_by_shard.items()}


def extract_records(
    candidate_groups: dict[tuple[str, str], CandidateGroup],
    shard_to_paths: dict[str, list[Path]],
) -> tuple[list[tuple[CandidateGroup, FastaRecord]], list[dict[str, str]]]:
    """Find selected candidate sequences in their source proteome paths."""

    target_ids_by_shard: dict[str, set[str]] = defaultdict(set)
    for proteome_shard, target_id in candidate_groups:
        target_ids_by_shard[proteome_shard].add(target_id)

    candidate_keys_by_path: dict[Path, set[tuple[str, str]]] = defaultdict(set)
    for proteome_shard, target_ids in target_ids_by_shard.items():
        paths = shard_to_paths[proteome_shard]
        path_by_accession: dict[str, Path] = {}
        for path in paths:
            accession = _proteome_accession_from_path(path)
            path_by_accession.setdefault(accession, path)
            path_by_accession.setdefault(_accession_without_version(accession), path)
        for target_id in target_ids:
            accession = _candidate_accession_from_target_id(target_id)
            exact_path = path_by_accession.get(accession) or path_by_accession.get(_accession_without_version(accession))
            if exact_path is not None:
                candidate_keys_by_path[exact_path].add((proteome_shard, target_id))
            else:
                for path in paths:
                    candidate_keys_by_path[path].add((proteome_shard, target_id))

    found: dict[tuple[str, str], FastaRecord] = {}
    missing_rows: list[dict[str, str]] = []
    for proteome_path in sorted(candidate_keys_by_path, key=lambda path: path.as_posix()):
        candidate_keys = candidate_keys_by_path[proteome_path]
        wanted_record_ids = {target_id for _, target_id in candidate_keys}
        wanted_keys_by_record_id: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for key in candidate_keys:
            wanted_keys_by_record_id[key[1]].append(key)
        for record in iter_fasta_records(proteome_path):
            if record.record_id not in wanted_record_ids:
                continue
            for key in wanted_keys_by_record_id[record.record_id]:
                if key in found:
                    raise ValueError(
                        f"P07 cannot disambiguate duplicate FASTA id {record.record_id!r} "
                        f"within proteome shard {key[0]!r}"
                    )
                found[key] = record

    for proteome_shard in sorted(target_ids_by_shard):
        paths = shard_to_paths[proteome_shard]
        for target_id in sorted(target_ids_by_shard[proteome_shard]):
            if (proteome_shard, target_id) in found:
                continue
            missing_rows.append(
                {
                    "proteome_shard": proteome_shard,
                    "target_id": target_id,
                    "expected_proteome_paths": _unique_join(path.as_posix() for path in paths),
                    "candidate_table_path": "",
                    "scan_manifest_path": "",
                    "notes": "Selected P06 candidate was absent from the proteome path(s) recorded in the P06 scan manifest.",
                }
            )

    ordered_found = [
        (candidate_groups[key], found[key])
        for key in sorted(found, key=lambda item: (item[0], item[1]))
    ]
    return ordered_found, missing_rows


def write_fasta_shards(
    found_records: list[tuple[CandidateGroup, FastaRecord]],
    fasta_dir: Path,
    *,
    candidate_table: Path,
    scan_manifest: Path,
    gtdb_release: str,
    generated_at_utc: str,
    sequences_per_shard: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Write deterministic P07 FASTA shards and their sequence manifest rows."""

    sequence_rows: list[dict[str, str]] = []
    shard_rows: list[dict[str, str]] = []
    used_sequence_ids: set[str] = set()
    for shard_index, chunk in enumerate(_chunks(found_records, sequences_per_shard), start=1):
        fasta_path = fasta_dir / f"p07_candidates_{shard_index:06d}.faa"
        with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
            for group, record in chunk:
                p07_sequence_id = _p07_sequence_id(group.proteome_shard, group.target_id, used_sequence_ids)
                expected_length = int(group.p06_target_length)
                if len(record.sequence) != expected_length:
                    raise ValueError(
                        f"P06 target_length mismatch for {group.proteome_shard}/{group.target_id}: "
                        f"P06={expected_length}, FASTA={len(record.sequence)}"
                    )
                handle.write(
                    f">{p07_sequence_id} original_id={group.target_id} "
                    f"proteome_shard={group.proteome_shard} "
                    f"families={_unique_join(sorted(group.families))}\n"
                )
                handle.write(_wrap_sequence(record.sequence))
                sequence_rows.append(
                    {
                        "p07_sequence_id": p07_sequence_id,
                        "proteome_shard": group.proteome_shard,
                        "target_id": group.target_id,
                        "source_proteome_path": record.source_path.as_posix(),
                        "target_length_from_p06": group.p06_target_length,
                        "sequence_length": str(len(record.sequence)),
                        "family_categories": _unique_join(sorted(group.families)),
                        "tiers": _unique_join(sorted(group.tiers)),
                        "max_full_sequence_score": f"{group.max_score:.1f}",
                        "p06_candidate_rows": str(group.p06_candidate_rows),
                        "fasta_shard": fasta_path.as_posix(),
                        "candidate_table_path": candidate_table.as_posix(),
                        "scan_manifest_path": scan_manifest.as_posix(),
                        "gtdb_release": gtdb_release,
                        "generated_at_utc": generated_at_utc,
                        "notes": (
                            "Sequence extracted from P03 GTDB proteome by P06 candidate id; "
                            "P07 annotation is domain/localization evidence only, not phenotype proof."
                        ),
                    }
                )
        shard_rows.append(
            {
                "fasta_shard": fasta_path.as_posix(),
                "sequence_count": str(len(chunk)),
            }
        )
    return sequence_rows, shard_rows


def build_command_manifest_rows(
    shard_rows: list[dict[str, str]],
    outdir: Path,
    *,
    interproscan_exe: str,
    interproscan_cpu: int,
    signalp_exe: str,
    signalp_mode: str,
    phobius_exe: str,
) -> list[dict[str, str]]:
    """Create planned-not-run annotation commands for each FASTA shard."""

    rows: list[dict[str, str]] = []
    for shard_row in shard_rows:
        fasta_path = Path(shard_row["fasta_shard"])
        shard_id = fasta_path.stem
        interpro_base = outdir / INTERPRO_DIRNAME / shard_id / "interproscan"
        rows.append(
            {
                "tool": "InterProScan",
                "fasta_shard": shard_id,
                "sequence_count": shard_row["sequence_count"],
                "input_fasta": fasta_path.as_posix(),
                "output_path": interpro_base.as_posix(),
                "command": (
                    f"{_shell_quote(interproscan_exe)} -i {_shell_quote(fasta_path.as_posix())} "
                    f"-f TSV,JSON,GFF3 -goterms -pa -cpu {interproscan_cpu} "
                    f"-b {_shell_quote(interpro_base.as_posix())}"
                ),
                "command_status": COMMAND_STATUS,
                "notes": (
                    "InterProScan scans protein-family/domain/site signatures; output supports architecture review "
                    "but does not prove degradation phenotype."
                ),
            }
        )

        signalp_dir = outdir / LOCALIZATION_DIRNAME / "signalp6" / shard_id
        rows.append(
            {
                "tool": "SignalP6",
                "fasta_shard": shard_id,
                "sequence_count": shard_row["sequence_count"],
                "input_fasta": fasta_path.as_posix(),
                "output_path": signalp_dir.as_posix(),
                "command": (
                    f"{_shell_quote(signalp_exe)} --fastafile {_shell_quote(fasta_path.as_posix())} "
                    f"--organism other --output_dir {_shell_quote(signalp_dir.as_posix())} "
                    f"--format none --mode {_shell_quote(signalp_mode)}"
                ),
                "command_status": COMMAND_STATUS,
                "notes": (
                    "SignalP6 'other' mode covers Bacteria and Archaea; signal-peptide prediction is localization "
                    "support, not functional phenotype evidence."
                ),
            }
        )

        if phobius_exe:
            phobius_path = outdir / LOCALIZATION_DIRNAME / "phobius" / f"{shard_id}.short.txt"
            rows.append(
                {
                    "tool": "Phobius",
                    "fasta_shard": shard_id,
                    "sequence_count": shard_row["sequence_count"],
                    "input_fasta": fasta_path.as_posix(),
                    "output_path": phobius_path.as_posix(),
                    "command": (
                        f"{_shell_quote(phobius_exe)} -short {_shell_quote(fasta_path.as_posix())} "
                        f"> {_shell_quote(phobius_path.as_posix())}"
                    ),
                    "command_status": COMMAND_STATUS,
                    "notes": (
                        "Optional local Phobius CLI must be verified on T141 before execution; it is used only for "
                        "signal-peptide/transmembrane topology review."
                    ),
                }
            )
    return rows


def summarize_p07_plan(
    candidate_groups: dict[tuple[str, str], CandidateGroup],
    sequence_rows: list[dict[str, str]],
    command_rows: list[dict[str, str]],
    *,
    selected_tiers: tuple[str, ...],
    missing_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Summarize the planned P07 annotation input set."""

    summary = [
        {"kind": "total", "name": "selected_unique_candidate_sequences", "count": str(len(candidate_groups))},
        {"kind": "total", "name": "extracted_sequences", "count": str(len(sequence_rows))},
        {"kind": "total", "name": "missing_sequences", "count": str(len(missing_rows))},
        {"kind": "total", "name": "fasta_shards", "count": str(len({row["fasta_shard"] for row in sequence_rows}))},
        {"kind": "total", "name": "command_rows", "count": str(len(command_rows))},
    ]
    for tier in selected_tiers:
        count = sum(1 for group in candidate_groups.values() if tier in group.tiers)
        summary.append({"kind": "selected_tier", "name": tier, "count": str(count)})
    for family, count in sorted(Counter(_iter_group_families(candidate_groups.values())).items()):
        summary.append({"kind": "family_category", "name": family, "count": str(count)})
    for tool, count in sorted(Counter(row["tool"] for row in command_rows).items()):
        summary.append({"kind": "planned_tool", "name": tool, "count": str(count)})
    return summary


def iter_fasta_records(path: Path) -> Iterator[FastaRecord]:
    """Yield FASTA records from plain or gzipped protein FASTA files."""

    with _open_text(path) as handle:
        header = ""
        sequence_parts: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield _fasta_record(header, sequence_parts, path)
                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if header:
            yield _fasta_record(header, sequence_parts, path)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _ensure_shards_have_paths(
    candidate_groups: dict[tuple[str, str], CandidateGroup],
    shard_to_paths: dict[str, list[Path]],
    scan_manifest: Path,
) -> None:
    missing = sorted({proteome_shard for proteome_shard, _ in candidate_groups if proteome_shard not in shard_to_paths})
    if missing:
        raise ValueError(
            f"{scan_manifest} does not contain proteome paths for selected P06 shard(s): {', '.join(missing[:10])}"
        )


def _candidate_accession_from_target_id(target_id: str) -> str:
    return target_id.split("|", 1)[0]


def _proteome_accession_from_path(path: Path) -> str:
    name = path.name
    for suffix in (".faa.gz", ".fa.gz", ".fasta.gz", ".faa", ".fa", ".fasta"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _accession_without_version(accession: str) -> str:
    prefix, separator, suffix = accession.rpartition(".")
    if separator and suffix.isdigit():
        return prefix
    return accession


def _require_fields(path: Path, fieldnames: list[str] | None, required: tuple[str, ...]) -> None:
    if fieldnames is None:
        raise ValueError(f"{path} is missing a tabular header")
    missing = [field_name for field_name in required if field_name not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _fasta_record(header: str, sequence_parts: list[str], source_path: Path) -> FastaRecord:
    record_id = header.split()[0]
    sequence = "".join(sequence_parts).replace(" ", "").upper()
    if not record_id:
        raise ValueError(f"{source_path} contains a FASTA record with an empty identifier")
    if not sequence:
        raise ValueError(f"{source_path} contains FASTA record {record_id!r} with an empty sequence")
    return FastaRecord(record_id=record_id, header=header, sequence=sequence, source_path=source_path)


def _chunks(items: list[tuple[CandidateGroup, FastaRecord]], size: int) -> Iterator[list[tuple[CandidateGroup, FastaRecord]]]:
    for index in range(0, len(items), size):
        yield items[index : index + size]


def _p07_sequence_id(proteome_shard: str, target_id: str, used_sequence_ids: set[str]) -> str:
    base = _safe_identifier(f"{proteome_shard}__{target_id}")
    sequence_id = base
    if sequence_id in used_sequence_ids:
        digest = hashlib.sha1(f"{proteome_shard}\0{target_id}".encode("utf-8")).hexdigest()[:10]
        sequence_id = f"{base}_{digest}"
    used_sequence_ids.add(sequence_id)
    return sequence_id


def _safe_identifier(value: str) -> str:
    allowed = []
    for char in value.strip():
        if char.isalnum() or char in {".", "_", "-", ":"}:
            allowed.append(char)
        else:
            allowed.append("_")
    return "".join(allowed) or "unnamed"


def _wrap_sequence(sequence: str, width: int = 60) -> str:
    return "".join(sequence[index : index + width] + "\n" for index in range(0, len(sequence), width))


def _unique_join(values: Iterable[str]) -> str:
    return ";".join(dict.fromkeys(values))


def _iter_group_families(groups: Iterable[CandidateGroup]) -> Iterator[str]:
    for group in groups:
        yield from group.families


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare P07 InterPro/localization annotation inputs from P06 candidates.")
    parser.add_argument("--candidate-table", type=Path, default=Path("05_hmmer_scan/p06_hmmer_candidates.tsv"))
    parser.add_argument("--scan-manifest", type=Path, default=Path("05_hmmer_scan/p06_hmmer_scan_manifest.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("06_domain_annotation"))
    parser.add_argument(
        "--include-tier",
        action="append",
        dest="include_tiers",
        default=None,
        help="P06 tier to include; repeat for multiple tiers. Defaults to High-confidence only.",
    )
    parser.add_argument("--sequences-per-shard", type=int, default=500)
    parser.add_argument("--gtdb-release", default="GTDB Release 11 R232")
    parser.add_argument("--interproscan-exe", default="interproscan.sh")
    parser.add_argument("--interproscan-cpu", type=int, default=8)
    parser.add_argument("--signalp-exe", default="signalp6")
    parser.add_argument("--signalp-mode", default="fast")
    parser.add_argument("--phobius-exe", default="")
    parser.add_argument("--allow-missing", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    outputs = prepare_p07_inputs(
        args.candidate_table,
        args.scan_manifest,
        args.outdir,
        include_tiers=args.include_tiers or DEFAULT_INCLUDE_TIERS,
        sequences_per_shard=args.sequences_per_shard,
        gtdb_release=args.gtdb_release,
        interproscan_exe=args.interproscan_exe,
        interproscan_cpu=args.interproscan_cpu,
        signalp_exe=args.signalp_exe,
        signalp_mode=args.signalp_mode,
        phobius_exe=args.phobius_exe,
        allow_missing=args.allow_missing,
    )
    print(f"Sequence manifest written: {outputs['sequence_manifest']}")
    print(f"Command manifest written: {outputs['command_manifest']}")
    print(f"Summary written: {outputs['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
