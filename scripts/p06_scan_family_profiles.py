"""Plan and parse P06 HMMER scans over GTDB proteome shards.

P06 keeps the raw HMMER `domtblout` outputs, then converts them into a
candidate catalog with conservative tiering. The scan step itself stays
explicitly sequence-based: it identifies homologous hits, but it does not
claim phenotype.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


OUTPUT_SCAN_MANIFEST_FILENAME = "p06_hmmer_scan_manifest.tsv"
OUTPUT_SCAN_SUMMARY_FILENAME = "p06_hmmer_scan_summary.tsv"
OUTPUT_CANDIDATE_TABLE_FILENAME = "p06_hmmer_candidates.tsv"
OUTPUT_CANDIDATE_SUMMARY_FILENAME = "p06_hmmer_candidate_summary.tsv"
RAW_DTBLOUT_DIRNAME = "raw_domtblout"
RAW_LOG_DIRNAME = "hmmer_logs"
COMMAND_STATUS = "planned_not_run"
MODEL_REGISTRY_REQUIRED_FIELDS = (
    "family_category",
    "approved_for_p06",
    "scan_permission",
    "model_path",
    "model_sha256",
    "model_specific_thresholds",
)
PROTEOME_FASTA_SUFFIXES = (
    ".faa.gz",
    ".fa.gz",
    ".fasta.gz",
    ".faa",
    ".fa",
    ".fasta",
)

SCAN_MANIFEST_FIELDNAMES = (
    "family_category",
    "proteome_shard",
    "proteome_count",
    "hmm_path",
    "model_sha256",
    "calibrated_full_score_threshold",
    "calibrated_hmm_coverage_threshold",
    "proteome_path",
    "domtblout_path",
    "main_output_path",
    "command",
    "command_status",
    "notes",
)

CANDIDATE_FIELDNAMES = (
    "family_category",
    "proteome_shard",
    "target_id",
    "target_accession",
    "query_id",
    "query_accession",
    "target_length",
    "query_length",
    "full_sequence_evalue",
    "full_sequence_score",
    "full_sequence_bias",
    "domain_index",
    "domain_count",
    "domain_c_evalue",
    "domain_i_evalue",
    "domain_score",
    "domain_bias",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "hmm_coverage",
    "target_coverage",
    "domain_overlap_fraction",
    "tier",
    "tier_reason",
    "domtblout_path",
    "calibrated_full_score_threshold",
    "calibrated_hmm_coverage_threshold",
)

SUMMARY_FIELDNAMES = ("kind", "name", "count")


def build_scan_manifest(
    hmm_dir: Path,
    proteome_dir: Path,
    outdir: Path,
    *,
    model_registry: Path,
    cpu: int = 1,
    proteomes_per_job: int = 1,
) -> dict[str, Path]:
    """Write a deterministic hmmsearch manifest for checksum-locked profiles only."""

    if proteomes_per_job < 1:
        raise ValueError("proteomes_per_job must be at least 1")

    approved_models = load_approved_hmm_models(model_registry, hmm_dir)
    proteome_paths = discover_proteome_paths(proteome_dir)
    proteome_chunks = _chunk_proteome_paths(proteome_paths, proteomes_per_job)
    raw_domtblout_dir = outdir / RAW_DTBLOUT_DIRNAME
    raw_log_dir = outdir / RAW_LOG_DIRNAME

    raw_domtblout_dir.mkdir(parents=True, exist_ok=True)
    raw_log_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str]] = []
    for hmm_path, model_sha256, score_threshold, hmm_coverage_threshold in approved_models:
        family_category = hmm_path.stem
        family_raw_dir = raw_domtblout_dir / _safe_identifier(family_category)
        family_log_dir = raw_log_dir / _safe_identifier(family_category)
        family_raw_dir.mkdir(parents=True, exist_ok=True)
        family_log_dir.mkdir(parents=True, exist_ok=True)

        for chunk_index, proteome_chunk in enumerate(proteome_chunks, start=1):
            shard_id = (
                _proteome_shard_id(proteome_chunk[0])
                if proteomes_per_job == 1
                else f"chunk_{chunk_index:06d}"
            )
            domtblout_path = family_raw_dir / f"{_safe_identifier(shard_id)}.domtblout"
            main_output_path = family_log_dir / f"{_safe_identifier(shard_id)}.txt"
            manifest_rows.append(
                {
                    "family_category": family_category,
                    "proteome_shard": shard_id,
                    "proteome_count": str(len(proteome_chunk)),
                    "hmm_path": _posix_path(hmm_path),
                    "model_sha256": model_sha256,
                    "calibrated_full_score_threshold": f"{score_threshold:.1f}",
                    "calibrated_hmm_coverage_threshold": f"{hmm_coverage_threshold:.6f}",
                    "proteome_path": _unique_join(_posix_path(path) for path in proteome_chunk),
                    "domtblout_path": _posix_path(domtblout_path),
                    "main_output_path": _posix_path(main_output_path),
                    "command": _hmmsearch_command(
                        hmm_path,
                        proteome_chunk,
                        domtblout_path,
                        main_output_path,
                        cpu=cpu,
                    ),
                    "command_status": COMMAND_STATUS,
                    "notes": (
                        "Prepared from the checksum-locked approved P05 HMM registry; hmmsearch was not run. "
                        "High-confidence parsing requires this model's calibrated score and HMM-coverage thresholds. "
                        "Raw domtblout and logs stay separate from derived candidate tables."
                    ),
                }
            )

    manifest_rows.sort(key=lambda row: (row["family_category"], row["proteome_shard"]))
    summary_rows = summarize_scan_manifest(manifest_rows)

    manifest_path = outdir / OUTPUT_SCAN_MANIFEST_FILENAME
    summary_path = outdir / OUTPUT_SCAN_SUMMARY_FILENAME
    write_tsv(manifest_path, manifest_rows, SCAN_MANIFEST_FIELDNAMES)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)
    return {"manifest": manifest_path, "summary": summary_path}


def summarize_scan_manifest(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize planned hmmsearch jobs."""

    family_count = len({row["family_category"] for row in rows})
    shard_count = len({row["proteome_shard"] for row in rows})
    summary = [
        {"kind": "total", "name": "scan_jobs", "count": str(len(rows))},
        {"kind": "total", "name": "families", "count": str(family_count)},
        {"kind": "total", "name": "proteome_shards", "count": str(shard_count)},
        {
            "kind": "command_status",
            "name": COMMAND_STATUS,
            "count": str(sum(row["command_status"] == COMMAND_STATUS for row in rows)),
        },
    ]
    return summary


def parse_scan_manifest(
    manifest_path: Path,
    outdir: Path,
) -> dict[str, Path]:
    """Convert completed domtblout files into a candidate catalog and summary."""

    scan_rows = load_scan_manifest(manifest_path)
    candidate_rows: list[dict[str, str]] = []
    missing_domtblout = 0
    for scan_row in scan_rows:
        domtblout_path = Path(scan_row["domtblout_path"])
        if not domtblout_path.is_file():
            missing_domtblout += 1
            continue
        raw_rows = load_domtblout_rows(domtblout_path)
        candidate_rows.extend(
            summarize_candidates_for_scan_row(scan_row, raw_rows, domtblout_path)
        )

    candidate_rows.sort(key=lambda row: (row["family_category"], row["proteome_shard"], row["target_id"]))
    summary_rows = summarize_candidates(candidate_rows, missing_domtblout=missing_domtblout)

    outdir.mkdir(parents=True, exist_ok=True)
    candidates_path = outdir / OUTPUT_CANDIDATE_TABLE_FILENAME
    summary_path = outdir / OUTPUT_CANDIDATE_SUMMARY_FILENAME
    write_tsv(candidates_path, candidate_rows, CANDIDATE_FIELDNAMES)
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDNAMES)
    return {"candidates": candidates_path, "summary": summary_path}


def summarize_candidates(
    candidate_rows: list[dict[str, str]],
    *,
    missing_domtblout: int = 0,
) -> list[dict[str, str]]:
    """Summarize the candidate catalog and tier balance."""

    summary = [
        {"kind": "total", "name": "candidate_rows", "count": str(len(candidate_rows))},
        {"kind": "total", "name": "families", "count": str(len({row["family_category"] for row in candidate_rows}))},
        {"kind": "total", "name": "proteome_shards", "count": str(len({row["proteome_shard"] for row in candidate_rows}))},
        {"kind": "total", "name": "missing_domtblout", "count": str(missing_domtblout)},
    ]
    tier_counts = Counter(row["tier"] for row in candidate_rows)
    for tier in ("High-confidence", "Review", "Rejected"):
        summary.append({"kind": "tier", "name": tier, "count": str(tier_counts.get(tier, 0))})
    return summary


def summarize_candidates_for_scan_row(
    scan_row: dict[str, str],
    raw_rows: list[dict[str, str]],
    domtblout_path: Path,
) -> list[dict[str, str]]:
    """Convert raw domtblout rows from one scan job into candidate rows."""

    if not raw_rows:
        return []

    domain_overlap_by_target = _domain_overlap_by_target(raw_rows)
    candidate_rows: list[dict[str, str]] = []
    for raw_row in raw_rows:
        candidate = _candidate_row_from_raw(scan_row, raw_row, domtblout_path)
        overlap = domain_overlap_by_target.get(candidate["target_id"], 0.0)
        candidate["domain_overlap_fraction"] = _format_fraction(overlap)
        tier, reason = classify_candidate(candidate)
        candidate["tier"] = tier
        candidate["tier_reason"] = reason
        candidate_rows.append(candidate)
    return candidate_rows


def classify_candidate(candidate: dict[str, str]) -> tuple[str, str]:
    """Assign a conservative tier to one parsed HMMER hit."""

    evalue = float(candidate["full_sequence_evalue"])
    score = float(candidate["full_sequence_score"])
    bias = float(candidate["full_sequence_bias"])
    hmm_coverage = float(candidate["hmm_coverage"])
    target_coverage = float(candidate["target_coverage"])
    domain_overlap = float(candidate["domain_overlap_fraction"])
    calibrated_score = float(candidate["calibrated_full_score_threshold"])
    calibrated_hmm_coverage = float(candidate["calibrated_hmm_coverage_threshold"])

    if (
        evalue <= 1e-5
        and score >= calibrated_score
        and hmm_coverage >= calibrated_hmm_coverage
        and target_coverage >= 0.20
        and bias <= 3.0
        and domain_overlap <= 0.35
    ):
        return "High-confidence", "passes this model's calibrated score and HMM-coverage thresholds plus review gates"

    if (
        evalue <= 1e-5
        and score >= 25.0
        and hmm_coverage >= 0.35
        and target_coverage >= 0.20
        and bias <= 3.0
        and domain_overlap <= 0.35
    ):
        return "Review", "passes relaxed score and coverage gates but needs manual review"

    return "Rejected", _rejection_reason(candidate)


def load_scan_manifest(path: Path) -> list[dict[str, str]]:
    """Load a scan manifest written by build_scan_manifest."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing_columns = [field for field in SCAN_MANIFEST_FIELDNAMES if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing_columns)}")

        rows: list[dict[str, str]] = []
        seen_pairs: set[tuple[str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            normalized.setdefault("proteome_count", "1")
            for field in SCAN_MANIFEST_FIELDNAMES:
                if not normalized.get(field):
                    raise ValueError(f"{path}:{line_number} missing required value for {field}")
            pair = (normalized["family_category"], normalized["proteome_shard"])
            if pair in seen_pairs:
                raise ValueError(
                    f"{path}:{line_number} has duplicate family_category/proteome_shard combination {pair!r}"
                )
            seen_pairs.add(pair)
            rows.append(normalized)
    return rows


def load_approved_hmm_models(model_registry: Path, hmm_dir: Path) -> list[tuple[Path, str, float, float]]:
    """Return only registry-approved HMMs after exact path and checksum checks."""

    if not model_registry.is_file():
        raise FileNotFoundError(f"P05 model registry not found: {model_registry}")
    if not hmm_dir.is_dir():
        raise FileNotFoundError(f"HMM directory not found: {hmm_dir}")

    with model_registry.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{model_registry} is missing a tabular header")
        missing_columns = [field for field in MODEL_REGISTRY_REQUIRED_FIELDS if field not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"{model_registry} is missing required columns: {', '.join(missing_columns)}")
        registry_rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]

    approved_rows = [
        row
        for row in registry_rows
        if row["approved_for_p06"].lower() == "yes" and row["scan_permission"].lower() == "approved"
    ]
    if not approved_rows:
        raise ValueError(
            f"{model_registry} has no P06-approved HMMs. Confirm the revised seed set, rebuild or retain models, "
            "record SHA256 values, and complete the calibration decision before planning a scan."
        )

    models: list[tuple[Path, str, float, float]] = []
    seen_families: set[str] = set()
    for line_number, row in enumerate(approved_rows, start=2):
        family = row["family_category"]
        expected_hash = row["model_sha256"].lower()
        model_path = Path(row["model_path"])
        if family in seen_families:
            raise ValueError(f"{model_registry}:{line_number} has duplicate approved family {family!r}")
        if not model_path.is_file():
            raise FileNotFoundError(f"{model_registry}:{line_number} approved HMM is missing: {model_path}")
        if model_path.parent.resolve() != hmm_dir.resolve():
            raise ValueError(f"{model_registry}:{line_number} approved HMM is outside --hmm-dir: {model_path}")
        observed_hash = _sha256(model_path)
        if observed_hash != expected_hash:
            raise ValueError(
                f"{model_registry}:{line_number} checksum mismatch for {model_path}: "
                f"expected {expected_hash}, observed {observed_hash}"
            )
        score_threshold, hmm_coverage_threshold = _parse_model_thresholds(
            row["model_specific_thresholds"], model_registry, line_number
        )
        seen_families.add(family)
        models.append((model_path, expected_hash, score_threshold, hmm_coverage_threshold))
    return sorted(models, key=lambda item: item[0].as_posix())


def _parse_model_thresholds(value: str, registry_path: Path, line_number: int) -> tuple[float, float]:
    """Read the two P05-derived acceptance thresholds from one approved registry row."""

    values: dict[str, float] = {}
    for token in value.split(";"):
        key, separator, raw_value = token.partition(">=")
        if not separator:
            continue
        try:
            parsed = float(raw_value)
        except ValueError as exc:
            raise ValueError(f"{registry_path}:{line_number} has an invalid calibrated threshold {token!r}") from exc
        if not math.isfinite(parsed) or parsed <= 0:
            raise ValueError(f"{registry_path}:{line_number} has a non-positive calibrated threshold {token!r}")
        values[key.strip()] = parsed
    try:
        return values["full_score"], values["hmm_coverage"]
    except KeyError as exc:
        raise ValueError(
            f"{registry_path}:{line_number} approved model must record full_score>=...;hmm_coverage>=..."
        ) from exc


def load_domtblout_rows(path: Path) -> list[dict[str, str]]:
    """Parse one HMMER domtblout file into rows."""

    rows: list[dict[str, str]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=22)
        if len(fields) < 23:
            raise ValueError(f"{path} has a malformed domtblout line: {line}")
        rows.append(
            {
                "target_id": fields[0],
                "target_accession": fields[1],
                "target_length": fields[2],
                "query_id": fields[3],
                "query_accession": fields[4],
                "query_length": fields[5],
                "full_sequence_evalue": fields[6],
                "full_sequence_score": fields[7],
                "full_sequence_bias": fields[8],
                "domain_index": fields[9],
                "domain_count": fields[10],
                "domain_c_evalue": fields[11],
                "domain_i_evalue": fields[12],
                "domain_score": fields[13],
                "domain_bias": fields[14],
                "hmm_from": fields[15],
                "hmm_to": fields[16],
                "ali_from": fields[17],
                "ali_to": fields[18],
                "env_from": fields[19],
                "env_to": fields[20],
                "accuracy": fields[21],
                "description": fields[22],
            }
        )
    return rows


def _candidate_row_from_raw(
    scan_row: dict[str, str],
    raw_row: dict[str, str],
    domtblout_path: Path,
) -> dict[str, str]:
    target_length = int(raw_row["target_length"])
    query_length = int(raw_row["query_length"])
    hmm_from = int(raw_row["hmm_from"])
    hmm_to = int(raw_row["hmm_to"])
    ali_from = int(raw_row["ali_from"])
    ali_to = int(raw_row["ali_to"])

    return {
        "family_category": scan_row["family_category"],
        "proteome_shard": scan_row["proteome_shard"],
        "target_id": raw_row["target_id"],
        "target_accession": raw_row["target_accession"],
        "query_id": raw_row["query_id"],
        "query_accession": raw_row["query_accession"],
        "target_length": raw_row["target_length"],
        "query_length": raw_row["query_length"],
        "full_sequence_evalue": raw_row["full_sequence_evalue"],
        "full_sequence_score": raw_row["full_sequence_score"],
        "full_sequence_bias": raw_row["full_sequence_bias"],
        "domain_index": raw_row["domain_index"],
        "domain_count": raw_row["domain_count"],
        "domain_c_evalue": raw_row["domain_c_evalue"],
        "domain_i_evalue": raw_row["domain_i_evalue"],
        "domain_score": raw_row["domain_score"],
        "domain_bias": raw_row["domain_bias"],
        "hmm_from": raw_row["hmm_from"],
        "hmm_to": raw_row["hmm_to"],
        "ali_from": raw_row["ali_from"],
        "ali_to": raw_row["ali_to"],
        "env_from": raw_row["env_from"],
        "env_to": raw_row["env_to"],
        "hmm_coverage": _format_fraction((hmm_to - hmm_from + 1) / query_length),
        "target_coverage": _format_fraction((ali_to - ali_from + 1) / target_length),
        "domain_overlap_fraction": "0.0000",
        "calibrated_full_score_threshold": scan_row["calibrated_full_score_threshold"],
        "calibrated_hmm_coverage_threshold": scan_row["calibrated_hmm_coverage_threshold"],
        "tier": "",
        "tier_reason": "",
        "domtblout_path": _posix_path(domtblout_path),
    }


def _domain_overlap_by_target(raw_rows: list[dict[str, str]]) -> dict[str, float]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in raw_rows:
        grouped[row["target_id"]].append((int(row["env_from"]), int(row["env_to"])))

    overlaps: dict[str, float] = {}
    for target_id, spans in grouped.items():
        if len(spans) < 2:
            overlaps[target_id] = 0.0
            continue
        max_overlap = 0.0
        for index, (start_a, end_a) in enumerate(spans):
            for start_b, end_b in spans[index + 1 :]:
                overlap = min(end_a, end_b) - max(start_a, start_b) + 1
                if overlap <= 0:
                    continue
                shorter = min(end_a - start_a + 1, end_b - start_b + 1)
                max_overlap = max(max_overlap, overlap / shorter)
        overlaps[target_id] = max_overlap
    return overlaps


def _rejection_reason(candidate: dict[str, str]) -> str:
    evalue = float(candidate["full_sequence_evalue"])
    score = float(candidate["full_sequence_score"])
    bias = float(candidate["full_sequence_bias"])
    hmm_coverage = float(candidate["hmm_coverage"])
    target_coverage = float(candidate["target_coverage"])
    domain_overlap = float(candidate["domain_overlap_fraction"])

    if evalue > 1e-5:
        return f"Rejected: full-sequence E-value {evalue:g} is above the relaxed scan gate"
    if score < 25.0:
        return f"Rejected: full-sequence bit score {score:g} is below the relaxed scan gate"
    if hmm_coverage < 0.35:
        return f"Rejected: HMM coverage {hmm_coverage:.3f} is below the relaxed scan gate"
    if target_coverage < 0.20:
        return f"Rejected: target coverage {target_coverage:.3f} is below the relaxed scan gate"
    if bias > 3.0:
        return f"Rejected: full-sequence bias {bias:g} is above the relaxed scan gate"
    if domain_overlap > 0.35:
        return f"Rejected: domain overlap fraction {domain_overlap:.3f} is above the relaxed scan gate"
    return "Rejected: hit did not meet the conservative review gate"


def discover_hmm_paths(hmm_dir: Path) -> list[Path]:
    """Discover custom HMMs in deterministic order."""

    if not hmm_dir.is_dir():
        raise FileNotFoundError(f"HMM directory not found: {hmm_dir}")
    return sorted(path for path in hmm_dir.iterdir() if path.is_file() and path.suffix == ".hmm")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_proteome_paths(proteome_dir: Path) -> list[Path]:
    """Discover GTDB proteome shards in deterministic order."""

    if not proteome_dir.is_dir():
        raise FileNotFoundError(f"Proteome directory not found: {proteome_dir}")
    return sorted(
        (
            path
            for path in proteome_dir.rglob("*")
            if path.is_file() and _has_proteome_fasta_suffix(path)
        ),
        key=lambda path: path.as_posix(),
    )


def _chunk_proteome_paths(paths: list[Path], proteomes_per_job: int) -> list[list[Path]]:
    return [paths[index : index + proteomes_per_job] for index in range(0, len(paths), proteomes_per_job)]


def _hmmsearch_command(
    hmm_path: Path,
    proteome_chunk: list[Path],
    domtblout_path: Path,
    main_output_path: Path,
    *,
    cpu: int,
) -> str:
    hmmsearch_part = (
        "hmmsearch --noali --acc --seed 42 --cpu "
        f"{cpu} --domtblout {_shell_quote(_posix_path(domtblout_path))} "
        f"-o {_shell_quote(_posix_path(main_output_path))} "
        f"{_shell_quote(_posix_path(hmm_path))}"
    )
    if len(proteome_chunk) == 1:
        return f"{hmmsearch_part} {_shell_quote(_posix_path(proteome_chunk[0]))}"
    zcat_inputs = " ".join(_shell_quote(_posix_path(path)) for path in proteome_chunk)
    return f"zcat {zcat_inputs} | {hmmsearch_part} -"


def _unique_join(values: Iterable[str]) -> str:
    return ";".join(values)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _format_fraction(value: float) -> str:
    if math.isnan(value) or math.isinf(value):
        raise ValueError("fraction must be finite")
    return f"{value:.4f}"


def _safe_identifier(value: str) -> str:
    allowed = []
    for char in value.strip():
        if char.isalnum() or char in {".", "_", "-", ":"}:
            allowed.append(char)
        else:
            allowed.append("_")
    result = "".join(allowed)
    return result or "unnamed"


def _has_proteome_fasta_suffix(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in PROTEOME_FASTA_SUFFIXES)


def _proteome_shard_id(path: Path) -> str:
    name = path.name
    lower_name = name.lower()
    for suffix in PROTEOME_FASTA_SUFFIXES:
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _posix_path(path: Path) -> str:
    return path.as_posix()


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plan and parse P06 HMMER scans over GTDB proteome shards.")
    parser.add_argument(
        "--hmm-dir",
        type=Path,
        default=Path("04_family_profiles/hmms"),
        help="Directory containing one custom HMM per family",
    )
    parser.add_argument(
        "--proteome-dir",
        type=Path,
        default=Path("03_gtdb_proteomes/faa"),
        help="Directory containing GTDB proteome shard FASTA files",
    )
    parser.add_argument(
        "--model-registry",
        type=Path,
        default=Path("04_family_profiles/manifests/p05_hmm_model_registry.tsv"),
        help="Checksum-locked P05 model registry; only rows approved for P06 are accepted",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("05_hmmer_scan"),
        help="Output directory for scan manifests and candidate tables",
    )
    parser.add_argument(
        "--proteomes-per-job",
        type=int,
        default=1,
        help="Number of proteome FASTA files to stream into each hmmsearch job",
    )
    parser.add_argument(
        "--parse-manifest",
        type=Path,
        default=Path("05_hmmer_scan/p06_hmmer_scan_manifest.tsv"),
        help="Existing scan manifest to parse into candidate tables",
    )
    parser.add_argument(
        "--parse-only",
        action="store_true",
        help="Parse an existing scan manifest instead of rebuilding it",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.parse_only:
        outputs = parse_scan_manifest(args.parse_manifest, args.outdir)
        print(f"Candidate table written: {outputs['candidates']}")
        print(f"Candidate summary written: {outputs['summary']}")
        return 0

    outputs = build_scan_manifest(
        args.hmm_dir,
        args.proteome_dir,
        args.outdir,
        model_registry=args.model_registry,
        proteomes_per_job=args.proteomes_per_job,
    )
    print(f"Scan manifest written: {outputs['manifest']}")
    print(f"Scan summary written: {outputs['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
