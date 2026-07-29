"""Prepare auditable P08 phylogeny inputs without running phylogeny tools.

This module joins accepted P06 candidates to P07 sequence/annotation records,
GTDB taxonomy, and checksum-verified P05 reference records.  It deliberately
does not create FASTA, command, alignment, or tree artifacts.
"""

from __future__ import annotations

import csv
import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.p02_select_benchmark_genomes import load_taxonomy, strip_gtdb_prefix


DEFAULT_INCLUDE_TIERS = ("High-confidence",)
REQUIRED_P07_TOOLS = ("InterProScan", "SignalP6")
EVIDENCE_BOUNDARY = "sequence_and_annotation_evidence_only_not_phenotype_proof"
BLOCK_FIELDS = ("stage", "family_category", "proteome_shard", "target_id", "reason", "source_path", "notes")
P06_FIELDS = ("family_category", "proteome_shard", "target_id", "target_accession", "target_length", "full_sequence_score", "hmm_coverage", "tier")
P07_SEQUENCE_FIELDS = ("p07_sequence_id", "proteome_shard", "target_id", "source_proteome_path", "target_length_from_p06", "sequence_length", "family_categories", "fasta_shard")
P07_STATUS_FIELDS = ("tool", "fasta_shard", "status")
REGISTRY_FIELDS = ("family_category", "approved_for_p06", "scan_permission", "model_sha256")
SEED_FIELDS = ("family_category", "seed_id", "source_accession", "sequence_path", "sequence_sha256")
CONTROL_FIELDS = ("family_category", "control_id", "control_role", "sequence_path", "sequence_sha256")


def is_approved_model(row: Mapping[str, str]) -> bool:
    """Return whether a P05 registry record is permitted for P06/P08 use."""
    return row["approved_for_p06"] == "yes" and row["scan_permission"] == "approved"


def route_family_size(sequence_count: int) -> str:
    """Return the P08 planning route for a candidate-family size."""
    if sequence_count < 200:
        return "mafft_linsi_then_review"
    if sequence_count <= 2000:
        return "mafft_auto_then_review"
    return "deterministic_representative_plan_then_fasttree_exploratory"


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _read_tsv(path: Path, required_fields: Sequence[str], label: str) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"{label} missing header")
            missing = [field for field in required_fields if field not in reader.fieldnames]
            if missing:
                raise ValueError(f"{label} missing required field: {', '.join(missing)}")
            rows = list(reader)
    except FileNotFoundError as error:
        raise ValueError(f"{label} source file not found: {path}") from error
    if not rows:
        raise ValueError(f"{label} empty table")
    for row_number, row in enumerate(rows, start=2):
        missing_values = [field for field in required_fields if not row.get(field, "").strip()]
        if missing_values:
            raise ValueError(f"{label} missing required value at row {row_number}: {', '.join(missing_values)}")
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _assembly_accession(source_proteome_path: str) -> str:
    name = Path(source_proteome_path).name
    for suffix in (".faa.gz", ".fa.gz", ".fasta.gz", ".faa", ".fa", ".fasta"):
        if name.endswith(suffix):
            return strip_gtdb_prefix(name[: -len(suffix)])
    return strip_gtdb_prefix(Path(name).stem)


def _blocked(blocks: list[dict[str, str]], reason: str, *, family_category: str = "", proteome_shard: str = "", target_id: str = "", source_path: str = "", notes: str = "") -> None:
    blocks.append({"stage": "P08", "family_category": family_category, "proteome_shard": proteome_shard, "target_id": target_id, "reason": reason, "source_path": source_path, "notes": notes})


def _fail(outdir: Path, blocks: list[dict[str, str]], message: str, **context: str) -> None:
    _blocked(blocks, message, **context)
    _write_tsv(outdir / "review" / "p08_blocked_records.tsv", BLOCK_FIELDS, sorted(blocks, key=lambda row: tuple(row[field] for field in BLOCK_FIELDS)))
    raise ValueError(message)


def prepare_p08_inputs(
    *,
    p06_candidate_table: Path,
    p07_sequence_table: Path,
    p07_status_table: Path,
    p05_model_registry: Path,
    p05_seed_table: Path,
    p05_control_table: Path,
    taxonomy_paths: Sequence[Path],
    outdir: Path,
    include_tiers: Sequence[str] = DEFAULT_INCLUDE_TIERS,
) -> dict[str, Path]:
    """Validate and write deterministic P08 manifests from existing P05--P07 data."""
    outdir = Path(outdir)
    blocks: list[dict[str, str]] = []
    include_tiers = tuple(include_tiers)
    if "Rejected" in include_tiers:
        _fail(outdir, blocks, "Rejected tier is not permitted by the P08 Python API")
    try:
        registry_rows = _read_tsv(Path(p05_model_registry), REGISTRY_FIELDS, "P05 registry")
        seed_rows = _read_tsv(Path(p05_seed_table), SEED_FIELDS, "P05 seed table")
        control_rows = _read_tsv(Path(p05_control_table), CONTROL_FIELDS, "P05 control table")
        p06_rows = _read_tsv(Path(p06_candidate_table), P06_FIELDS, "P06 candidate table")
        p07_rows = _read_tsv(Path(p07_sequence_table), P07_SEQUENCE_FIELDS, "P07 sequence table")
        status_rows = _read_tsv(Path(p07_status_table), P07_STATUS_FIELDS, "P07 status table")
    except ValueError as error:
        _fail(outdir, blocks, str(error))

    p06_keys: set[tuple[str, str, str]] = set()
    for row in p06_rows:
        key = (row["family_category"], row["proteome_shard"], row["target_id"])
        if key in p06_keys:
            _fail(outdir, blocks, "duplicate identical P06 row key", family_category=key[0], proteome_shard=key[1], target_id=key[2], source_path=str(p06_candidate_table))
        p06_keys.add(key)
    selected = [row for row in p06_rows if row["tier"] in include_tiers]
    if not selected:
        _fail(outdir, blocks, "P06 candidate table has no selected tiers", source_path=str(p06_candidate_table))

    approved_models = {row["family_category"]: row for row in registry_rows if is_approved_model(row)}
    registry_families = {row["family_category"] for row in registry_rows}
    for row in selected:
        family = row["family_category"]
        if family not in registry_families:
            _fail(outdir, blocks, "unknown P06 family", family_category=family, proteome_shard=row["proteome_shard"], target_id=row["target_id"], source_path=str(p06_candidate_table))
        if family not in approved_models:
            _fail(outdir, blocks, "missing approved P05 model", family_category=family, proteome_shard=row["proteome_shard"], target_id=row["target_id"], source_path=str(p05_model_registry))

    p07_by_key = {(row["proteome_shard"], row["target_id"]): row for row in p07_rows}
    statuses = {(row["tool"], row["fasta_shard"]): row["status"] for row in status_rows}
    try:
        taxonomy = load_taxonomy([Path(path) for path in taxonomy_paths])
    except (OSError, ValueError) as error:
        _fail(outdir, blocks, f"taxonomy loading failed: {error}", source_path=", ".join(map(str, taxonomy_paths)))

    references_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for record_kind, rows, identifier_field in (("seed", seed_rows, "seed_id"), ("control", control_rows, "control_id")):
        for row in rows:
            path = Path(row["sequence_path"])
            try:
                verified = _sha256(path)
            except OSError as error:
                _fail(outdir, blocks, f"sequence file unreadable: {path}", family_category=row["family_category"], source_path=str(path), notes=str(error))
            if verified != row["sequence_sha256"]:
                _fail(outdir, blocks, "SHA-256 mismatch", family_category=row["family_category"], source_path=str(path), notes=f"{identifier_field}={row[identifier_field]}")
            references_by_family[row["family_category"]].append({
                "family_category": row["family_category"], "record_kind": record_kind, "record_id": row[identifier_field], "source_accession": row.get("source_accession", ""), "control_role": row.get("control_role", ""), "sequence_path": str(path), "sequence_sha256": row["sequence_sha256"], "verified_sha256": verified, "evidence": row.get("evidence", ""), "notes": row.get("notes", ""),
            })

    candidate_rows: list[dict[str, str]] = []
    taxonomy_rows: list[dict[str, str]] = []
    for p06 in selected:
        family = p06["family_category"]
        family_references = references_by_family[family]
        if not any(row["record_kind"] == "seed" for row in family_references):
            _fail(outdir, blocks, "missing P05 seed", family_category=family, proteome_shard=p06["proteome_shard"], target_id=p06["target_id"])
        if not any(row["record_kind"] == "control" for row in family_references):
            _fail(outdir, blocks, "missing P05 control", family_category=family, proteome_shard=p06["proteome_shard"], target_id=p06["target_id"])
        key = (p06["proteome_shard"], p06["target_id"])
        p07 = p07_by_key.get(key)
        if p07 is None:
            _fail(outdir, blocks, "missing P07 sequence match", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table))
        if p06["target_length"] != p07["target_length_from_p06"] or p06["target_length"] != p07["sequence_length"]:
            _fail(outdir, blocks, "P06/P07 sequence length mismatch", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table))
        assembly_accession = _assembly_accession(p07["source_proteome_path"])
        lineage = taxonomy.get(assembly_accession)
        if not lineage:
            _fail(outdir, blocks, "taxonomy absence blocks processing", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["source_proteome_path"])
        shard_stem = Path(p07["fasta_shard"]).stem
        annotation_complete = all(statuses.get((tool, shard_stem)) in {"completed", "skipped_existing"} for tool in REQUIRED_P07_TOOLS)
        candidate = {
            **{field: p06[field] for field in P06_FIELDS},
            "p07_sequence_id": p07["p07_sequence_id"], "source_proteome_path": p07["source_proteome_path"], "fasta_shard": p07["fasta_shard"], "p07_annotation_status": "completed" if annotation_complete else "incomplete", "assembly_accession": assembly_accession, "taxonomy_lineage": lineage, "model_sha256": approved_models[family]["model_sha256"], "evidence_boundary": EVIDENCE_BOUNDARY,
        }
        candidate_rows.append(candidate)
        taxonomy_rows.append({"family_category": family, "proteome_shard": p06["proteome_shard"], "target_id": p06["target_id"], "assembly_accession": assembly_accession, "taxonomy_lineage": lineage})

    sort_key = lambda row: (row["family_category"], row["proteome_shard"], row["target_id"])
    candidate_rows.sort(key=sort_key)
    taxonomy_rows.sort(key=sort_key)
    reference_rows = sorted((row for family in {row["family_category"] for row in selected} for row in references_by_family[family]), key=lambda row: (row["family_category"], row["record_kind"], row["record_id"]))
    summary_rows = [{"family_category": family, "candidate_count": str(sum(row["family_category"] == family for row in candidate_rows)), "route": route_family_size(sum(row["family_category"] == family for row in candidate_rows))} for family in sorted({row["family_category"] for row in candidate_rows})]

    outputs = {
        "candidate_manifest": outdir / "manifests" / "p08_candidate_manifest.tsv",
        "taxonomy_join": outdir / "gtdb_mapping" / "p08_taxonomy_join.tsv",
        "family_reference_manifest": outdir / "manifests" / "p08_family_reference_manifest.tsv",
        "preparation_summary": outdir / "manifests" / "p08_preparation_summary.tsv",
    }
    _write_tsv(outputs["candidate_manifest"], tuple(candidate_rows[0]), candidate_rows)
    _write_tsv(outputs["taxonomy_join"], ("family_category", "proteome_shard", "target_id", "assembly_accession", "taxonomy_lineage"), taxonomy_rows)
    _write_tsv(outputs["family_reference_manifest"], ("family_category", "record_kind", "record_id", "source_accession", "control_role", "sequence_path", "sequence_sha256", "verified_sha256", "evidence", "notes"), reference_rows)
    _write_tsv(outputs["preparation_summary"], ("family_category", "candidate_count", "route"), summary_rows)
    return outputs
