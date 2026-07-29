"""Prepare auditable P08 phylogeny inputs without running phylogeny tools.

This module joins accepted P06 candidates to P07 sequence/annotation records,
GTDB taxonomy, and checksum-verified P05 reference records.  It deliberately
does not execute alignment, tree, or other external phylogeny tools.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.p02_select_benchmark_genomes import load_taxonomy, strip_gtdb_prefix


DEFAULT_INCLUDE_TIERS = ("High-confidence",)
REQUIRED_P07_TOOLS = ("InterProScan", "SignalP6")
EVIDENCE_BOUNDARY = "sequence_and_annotation_evidence_only_not_phenotype_proof"
BLOCK_FIELDS = ("stage", "family_category", "proteome_shard", "target_id", "reason", "source_path", "notes")
P06_FIELDS = (
    "family_category", "proteome_shard", "target_id", "target_accession", "target_length",
    "full_sequence_score", "hmm_coverage", "calibrated_full_score_threshold",
    "calibrated_hmm_coverage_threshold", "tier",
)
P07_SEQUENCE_FIELDS = (
    "p07_sequence_id", "proteome_shard", "target_id", "source_proteome_path",
    "target_length_from_p06", "sequence_length", "family_categories", "fasta_shard",
    "candidate_table_path", "scan_manifest_path", "gtdb_release",
)
P07_STATUS_FIELDS = ("tool", "fasta_shard", "input_fasta", "output_path", "status")
P03_MANIFEST_FIELDS = ("accession", "faa_path", "status")
P03_QC_FIELDS = ("accession", "faa_path", "status")
REGISTRY_FIELDS = ("family_category", "approved_for_p06", "scan_permission", "model_sha256")
SEED_FIELDS = ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path", "sequence_sha256")
CONTROL_FIELDS = ("family_category", "model_sha256", "control_id", "control_role", "sequence_path", "sequence_sha256", "source_checksum_kind")
CORE_FAMILY = "extracellular_pha_depolymerase_core"
CORE_SEED_REGISTRY_FILENAME = "p05_extracellular_core_seed_registry.tsv"
CORE_CLOSE_CONTROLS_FILENAME = "p05_extracellular_core_close_controls.tsv"
LEGACY_CALIBRATION_CONTROL_PANEL_FILENAME = "p05_hmm_calibration_control_panel.tsv"
CORE_SEED_FIELDS = ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path")
CORE_CLOSE_CONTROL_FIELDS = ("source_accession", "sequence_path", "notes")
FAMILY_INPUT_FIELDS = (
    "family_category", "record_kind", "record_identity", "input_fasta_path", "input_sha256",
    "source_path", "source_sha256", "source_checksum_kind", "model_sha256", "source_model_sha256",
    "model_provenance", "model_provenance_source_path", "p06_proteome_shard", "p06_target_id",
    "p07_sequence_id", "source_accession_or_assembly", "evidence_or_control_role",
    "is_gtdb_candidate", "evidence_boundary",
)
COMMAND_FIELDS = (
    "family_category", "command_status", "input_fasta_path", "input_sha256",
    "candidate_input_record_count", "total_input_record_count", "route", "alignment_fasta_path",
    "representative_input_fasta_path", "fasttree_tree_path", "iqtree_prefix", "representative_plan",
    "representative_selection_algorithm", "representative_selection_algorithm_version",
    "representative_selection_parameters", "representative_selection_mapping_path",
    "representative_selection_mapping_sha256", "representative_selection_mapping_record_count",
    "representative_materialization_status",
    "mafft_template", "fasttree_template", "iqtree2_template", "iqtree2_annotation",
    "rooting_policy", "evidence_boundary",
)
SUMMARY_FIELDS = (
    "family_category", "candidate_count", "route", "family_fasta_path", "family_fasta_sha256",
    "total_fasta_record_count",
)
TAXONOMY_JOIN_FIELDS = (
    "family_category", "proteome_shard", "target_id", "assembly_accession", "taxonomy_lineage",
    "taxonomy_source_role", "taxonomy_source_path", "taxonomy_source_sha256",
)
INPUT_PROVENANCE_FIELDS = ("input_role", "input_path", "input_sha256", "input_usage")
ROOTING_POLICY = "explicit_accessioned_outgroup_required; otherwise midpoint_display_only"
IQTREE2_TEMPLATE = "iqtree2 -s {alignment_fasta} -m TEST -B 1000 --prefix {iqtree_prefix}"
CORE_SEED_COUNT = 17
CORE_CROSS_FAMILY_CHALLENGE_COUNT = 15
CORE_CLOSE_CONTROL_COUNT = 5
CORE_HARD_PANEL_COUNT = CORE_CROSS_FAMILY_CHALLENGE_COUNT + CORE_CLOSE_CONTROL_COUNT
REPRESENTATIVE_SELECTION_ALGORITHM = "cluster_then_sha256_tiebreak"
REPRESENTATIVE_SELECTION_ALGORITHM_VERSION = "v1"
REPRESENTATIVE_SELECTION_PARAMETERS = "cluster_identity=0.99;representative_tiebreak=sequence_sha256_then_record_identity;input_order=record_identity;target_count=requires_separate_approval"


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


def planned_command_templates(
    candidate_count: int,
    *,
    mafft_exe: str = "mafft",
    iqtree_exe: str = "iqtree2",
    fasttree_exe: str = "FastTree",
) -> dict[str, str]:
    """Return unexecuted command templates appropriate to a candidate count."""
    route = route_family_size(candidate_count)
    templates = {
        "route": route,
        "mafft_template": "",
        "fasttree_template": "",
        "representative_plan": "",
        "representative_selection_algorithm": "",
        "representative_selection_algorithm_version": "",
        "representative_selection_parameters": "",
        "representative_selection_mapping_sha256": "",
        "representative_selection_mapping_record_count": "",
        "representative_materialization_status": "not_applicable",
        "iqtree2_template": IQTREE2_TEMPLATE.replace("iqtree2", iqtree_exe, 1),
    }
    if route == "mafft_linsi_then_review":
        templates["mafft_template"] = f"{mafft_exe} --localpair --maxiterate 1000 --thread {{threads}} --inputorder {{input_fasta}} > {{alignment_fasta}}"
    elif route == "mafft_auto_then_review":
        templates["mafft_template"] = f"{mafft_exe} --auto --thread {{threads}} --inputorder {{input_fasta}} > {{alignment_fasta}}"
    else:
        templates["representative_plan"] = "deterministic representative selection is declared but not materialized; a separately approved selection mapping is required before exploratory FastTree"
        templates["representative_selection_algorithm"] = REPRESENTATIVE_SELECTION_ALGORITHM
        templates["representative_selection_algorithm_version"] = REPRESENTATIVE_SELECTION_ALGORITHM_VERSION
        templates["representative_selection_parameters"] = REPRESENTATIVE_SELECTION_PARAMETERS
        templates["representative_selection_mapping_sha256"] = "not_materialized"
        templates["representative_selection_mapping_record_count"] = "0"
        templates["representative_materialization_status"] = "requires_separate_approval"
        templates["fasttree_template"] = f"{fasttree_exe} -lg {{representative_alignment_fasta}} > {{fasttree_tree}}"
    return templates


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


def _load_control_rows(path: Path) -> list[dict[str, str]]:
    """Load controls with explicit checksum kinds and a named legacy-panel adapter."""
    rows = _read_tsv(path, CONTROL_FIELDS[:-1], "P05 control table")
    for row_number, row in enumerate(rows, start=2):
        checksum_kind = row.get("source_checksum_kind", "").strip()
        if checksum_kind:
            if checksum_kind not in {"file_sha256", "residue_sha256"}:
                raise ValueError(f"P05 control table invalid source_checksum_kind at row {row_number}: {checksum_kind}")
            continue
        if path.name == LEGACY_CALIBRATION_CONTROL_PANEL_FILENAME:
            row["source_checksum_kind"] = "residue_sha256"
            continue
        raise ValueError(f"P05 control table missing required value at row {row_number}: source_checksum_kind")
    return rows


def _unique_p03_rows(rows: Sequence[dict[str, str]], label: str, outdir: Path, blocks: list[dict[str, str]], source_path: Path) -> dict[str, dict[str, str]]:
    """Index accepted P03 rows by accession without permitting silent replacement."""
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row["accession"]
        if accession in indexed:
            _fail(outdir, blocks, f"duplicate {label} accession", source_path=str(source_path), notes=f"accession={accession}")
        if row["status"] not in {"ok", "completed"}:
            _fail(outdir, blocks, f"{label} row is not accepted", source_path=str(source_path), notes=f"accession={accession}; status={row['status']}")
        indexed[accession] = row
    return indexed


def validate_tracked_core_authority_tables(repo_root: Path) -> dict[str, object]:
    """Validate the compact tracked core authority contract without reading raw sequences."""
    manifests = Path(repo_root) / "04_family_profiles" / "manifests"
    registry_path = manifests / "p05_hmm_model_registry.tsv"
    ordinary_seeds_path = manifests / "p05_hmm_seed_registry.tsv"
    ordinary_controls_path = manifests / LEGACY_CALIBRATION_CONTROL_PANEL_FILENAME
    seeds_path = manifests / CORE_SEED_REGISTRY_FILENAME
    controls_path = manifests / CORE_CLOSE_CONTROLS_FILENAME
    registry_rows = _read_tsv(registry_path, REGISTRY_FIELDS, "tracked P05 model registry")
    core_rows = [row for row in registry_rows if row["family_category"] == CORE_FAMILY]
    if len(core_rows) != 1 or not is_approved_model(core_rows[0]):
        raise ValueError("tracked core model registry must contain one approved core row")
    ordinary_seed_rows = _read_tsv(ordinary_seeds_path, SEED_FIELDS, "tracked P05 seed registry")
    ordinary_control_rows = _read_tsv(ordinary_controls_path, CONTROL_FIELDS[:-1], "tracked P05 control panel")
    if any(row["family_category"] == CORE_FAMILY for row in ordinary_seed_rows) or any(
        row["family_category"] == CORE_FAMILY for row in ordinary_control_rows
    ):
        raise ValueError("tracked ordinary P05 tables must not contain direct core rows")
    seed_rows = _read_tsv(seeds_path, CORE_SEED_FIELDS, "tracked extracellular core seed registry")
    close_rows = _read_tsv(controls_path, CORE_CLOSE_CONTROL_FIELDS, "tracked extracellular core close controls")
    if len(seed_rows) != CORE_SEED_COUNT or len({row["source_accession"] for row in seed_rows}) != CORE_SEED_COUNT:
        raise ValueError("tracked core seed authority must contain 17 unique accessions")
    if len(close_rows) != CORE_CLOSE_CONTROL_COUNT or len({row["source_accession"] for row in close_rows}) != CORE_CLOSE_CONTROL_COUNT:
        raise ValueError("tracked core close-control authority must contain 5 unique accessions")
    if any(row["family_category"] != CORE_FAMILY or row["model_sha256"] != core_rows[0]["model_sha256"] for row in seed_rows):
        raise ValueError("tracked core seed authority model contract mismatch")
    approved_models = {row["family_category"]: row for row in registry_rows if is_approved_model(row)}
    cross_family_rows = [
        row for row in ordinary_seed_rows
        if row["family_category"] in approved_models
        and row["family_category"] != CORE_FAMILY
        and row["model_sha256"] == approved_models[row["family_category"]]["model_sha256"]
    ]
    if (
        len(cross_family_rows) != CORE_CROSS_FAMILY_CHALLENGE_COUNT
        or len({row["source_accession"] for row in cross_family_rows}) != CORE_CROSS_FAMILY_CHALLENGE_COUNT
    ):
        raise ValueError("tracked core cross-family authority must contain 15 unique challenges")
    return {
        "core_seed_count": len(seed_rows),
        "cross_family_challenge_count": len(cross_family_rows),
        "close_control_count": len(close_rows),
        "hard_panel_count": len(cross_family_rows) + len(close_rows),
        "ordinary_seed_registry_path": ordinary_seeds_path,
        "ordinary_seed_registry_sha256": _sha256(ordinary_seeds_path),
        "ordinary_control_panel_path": ordinary_controls_path,
        "ordinary_control_panel_sha256": _sha256(ordinary_controls_path),
        "core_seed_registry_path": seeds_path,
        "core_seed_registry_sha256": _sha256(seeds_path),
        "close_controls_path": controls_path,
        "close_controls_sha256": _sha256(controls_path),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalized_path(path: str | Path) -> str:
    """Return a stable local path identity without requiring the target to exist."""
    raw_path = Path(path)
    if not str(raw_path).strip():
        raise ValueError("path is empty")
    return os.path.normcase(os.path.normpath(str(raw_path.resolve(strict=False))))


def _same_normalized_path(left: str | Path, right: str | Path) -> bool:
    return _normalized_path(left) == _normalized_path(right)


def _load_taxonomy_records(
    paths: Sequence[Path],
    source_roles: Sequence[str] | None,
) -> dict[str, dict[str, str]]:
    """Load GTDB taxonomy while retaining CLI source-role provenance."""
    paths = tuple(Path(path) for path in paths)
    if source_roles is None:
        combined = load_taxonomy(paths)
        return {
            accession: {
                "lineage": lineage,
                "source_role": "combined_taxonomy_input",
                "source_path": "",
                "source_sha256": "",
            }
            for accession, lineage in combined.items()
        }
    if len(paths) != len(source_roles):
        raise ValueError("taxonomy paths and source roles must have the same length")
    expected_domains = {
        "bac120_taxonomy": "d__Bacteria",
        "ar53_taxonomy": "d__Archaea",
    }
    records: dict[str, dict[str, str]] = {}
    for path, role in zip(paths, source_roles):
        expected_domain = expected_domains.get(role)
        if expected_domain is None:
            raise ValueError(f"unsupported taxonomy source role: {role}")
        source_sha256 = _sha256(path)
        for accession, lineage in load_taxonomy([path]).items():
            if lineage.split(";", maxsplit=1)[0] != expected_domain:
                label = "Bac120 taxonomy" if role == "bac120_taxonomy" else "Ar53 taxonomy"
                raise ValueError(f"{label} contains a lineage outside {expected_domain}: {accession}")
            previous = records.get(accession)
            if previous is not None and previous["source_role"] != role:
                raise ValueError(
                    f"taxonomy accession occurs in both sources: {accession}; "
                    f"{previous['source_role']} and {role}"
                )
            records[accession] = {
                "lineage": lineage,
                "source_role": role,
                "source_path": str(path),
                "source_sha256": source_sha256,
            }
    return records


def _read_fasta(path: Path) -> dict[str, str]:
    """Read a simple unaligned protein FASTA and reject ambiguous records."""
    records: dict[str, str] = {}
    current_id: str | None = None
    sequence_lines: list[str] = []

    def finish_record() -> None:
        if current_id is None:
            return
        sequence = "".join(sequence_lines)
        if not sequence or any(not (residue.isalpha() or residue == "*") for residue in sequence):
            raise ValueError("invalid or empty sequence")
        if current_id in records:
            raise ValueError(f"duplicate FASTA identifier: {current_id}")
        records[current_id] = sequence.upper()

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                current_id = header.split(maxsplit=1)[0] if header else ""
                if not current_id:
                    raise ValueError("empty FASTA identifier")
                sequence_lines = []
            elif not line.strip():
                continue
            elif current_id is None:
                raise ValueError("sequence data before FASTA header")
            else:
                sequence_lines.append(line.strip())
    finish_record()
    if not records:
        raise ValueError("no FASTA records")
    return records


def _residue_sha256(sequence: str) -> str:
    """Hash normalized residues for legacy P05 calibration-control records."""
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _fasta_safe(value: str) -> str:
    """Return a deterministic non-opaque FASTA-header component."""
    cleaned = "".join(character if character.isalnum() or character in "_.-" else "-" for character in value.strip())
    return cleaned or "unknown"


def _write_fasta(path: Path, records: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            header = "|".join((
                _fasta_safe(record["record_kind"]),
                _fasta_safe(record["record_id"]),
                _fasta_safe(record["family_category"]),
                _fasta_safe(record["source_accession_or_assembly"]),
            ))
            handle.write(f">{header}\n{record['sequence']}\n")


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
    mafft_exe: str = "mafft",
    iqtree_exe: str = "iqtree2",
    fasttree_exe: str = "FastTree",
    taxonomy_source_roles: Sequence[str] | None = None,
    additional_provenance_inputs: Mapping[str, Path] | None = None,
    gtdb_release: str | None = None,
    p06_scan_manifest: Path | None = None,
    p03_prediction_manifest: Path | None = None,
    p03_prediction_qc: Path | None = None,
) -> dict[str, Path]:
    """Validate and write deterministic P08 manifests from existing P05--P07 data."""
    outdir = Path(outdir)
    blocks: list[dict[str, str]] = []
    include_tiers = tuple(include_tiers)
    if "Rejected" in include_tiers:
        _fail(outdir, blocks, "Rejected tier is not permitted by the P08 Python API")
    provenance_paths = {
        "p06_scan_manifest": p06_scan_manifest,
        "p03_prediction_manifest": p03_prediction_manifest,
        "p03_prediction_qc": p03_prediction_qc,
    }
    missing_provenance = [role for role, path in provenance_paths.items() if path is None]
    if not gtdb_release or not gtdb_release.strip() or missing_provenance:
        _fail(outdir, blocks, "P08 provenance chain requires explicit GTDB release, P06 scan manifest, P03 prediction manifest, and P03 prediction QC")
    provenance_paths = {role: Path(path) for role, path in provenance_paths.items()}
    for role, path in provenance_paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            _fail(outdir, blocks, f"P08 provenance chain requires an existing nonempty {role}", source_path=str(path))
    try:
        registry_rows = _read_tsv(Path(p05_model_registry), REGISTRY_FIELDS, "P05 registry")
        seed_rows = _read_tsv(Path(p05_seed_table), SEED_FIELDS, "P05 seed table")
        control_rows = _load_control_rows(Path(p05_control_table))
        p06_rows = _read_tsv(Path(p06_candidate_table), P06_FIELDS, "P06 candidate table")
        p07_rows = _read_tsv(Path(p07_sequence_table), P07_SEQUENCE_FIELDS, "P07 sequence table")
        status_rows = _read_tsv(Path(p07_status_table), P07_STATUS_FIELDS, "P07 status table")
        p03_manifest_rows = _read_tsv(provenance_paths["p03_prediction_manifest"], P03_MANIFEST_FIELDS, "P03 prediction manifest")
        p03_qc_rows = _read_tsv(provenance_paths["p03_prediction_qc"], P03_QC_FIELDS, "P03 prediction QC")
    except ValueError as error:
        _fail(outdir, blocks, str(error))
    p03_manifest_by_accession = _unique_p03_rows(
        p03_manifest_rows, "P03 prediction manifest", outdir, blocks, provenance_paths["p03_prediction_manifest"]
    )
    p03_qc_by_accession = _unique_p03_rows(
        p03_qc_rows, "P03 prediction QC", outdir, blocks, provenance_paths["p03_prediction_qc"]
    )

    p06_keys: set[tuple[str, str, str]] = set()
    for row in p06_rows:
        key = (row["family_category"], row["proteome_shard"], row["target_id"])
        if key in p06_keys:
            _fail(outdir, blocks, "duplicate identical P06 row key", family_category=key[0], proteome_shard=key[1], target_id=key[2], source_path=str(p06_candidate_table))
        p06_keys.add(key)
    approved_models = {row["family_category"]: row for row in registry_rows if is_approved_model(row)}
    registry_families = {row["family_category"] for row in registry_rows}
    for row in p06_rows:
        family = row["family_category"]
        if family not in registry_families:
            _fail(outdir, blocks, "unknown P06 family", family_category=family, proteome_shard=row["proteome_shard"], target_id=row["target_id"], source_path=str(p06_candidate_table))
        if family not in approved_models:
            _fail(outdir, blocks, "missing approved P05 model", family_category=family, proteome_shard=row["proteome_shard"], target_id=row["target_id"], source_path=str(p05_model_registry))
    selected = [row for row in p06_rows if row["tier"] in include_tiers]
    if not selected:
        _fail(outdir, blocks, "P06 candidate table has no selected tiers", source_path=str(p06_candidate_table))
    selected_families = {row["family_category"] for row in selected}
    if CORE_FAMILY in selected_families:
        if any(row["family_category"] == CORE_FAMILY for row in seed_rows):
            _fail(
                outdir,
                blocks,
                "direct P05 core seed rows are prohibited; derive the core from its authoritative registry",
                family_category=CORE_FAMILY,
                source_path=str(p05_seed_table),
            )
        if any(row["family_category"] == CORE_FAMILY for row in control_rows):
            _fail(
                outdir,
                blocks,
                "direct P05 core control rows are prohibited; derive the hard panel from its authoritative sources",
                family_category=CORE_FAMILY,
                source_path=str(p05_control_table),
            )

    expected_gtdb_release = gtdb_release.strip()
    expected_candidate_table = Path(p06_candidate_table)
    expected_scan_manifest = provenance_paths["p06_scan_manifest"]
    expected_candidate_table_sha256 = _sha256(expected_candidate_table)
    expected_scan_manifest_sha256 = _sha256(expected_scan_manifest)

    p07_by_key: dict[tuple[str, str], dict[str, str]] = {}
    p07_ids: set[str] = set()
    for row in p07_rows:
        key = (row["proteome_shard"], row["target_id"])
        if row["gtdb_release"].strip() != expected_gtdb_release:
            _fail(
                outdir,
                blocks,
                "P07 GTDB release mismatch",
                proteome_shard=key[0],
                target_id=key[1],
                source_path=str(p07_sequence_table),
                notes=f"p07_gtdb_release={row['gtdb_release']}; explicit_gtdb_release={expected_gtdb_release}",
            )
        if not _same_normalized_path(row["candidate_table_path"], expected_candidate_table):
            _fail(
                outdir,
                blocks,
                "P07 candidate_table_path mismatch",
                proteome_shard=key[0],
                target_id=key[1],
                source_path=str(p07_sequence_table),
                notes=f"p07_candidate_table_path={row['candidate_table_path']}; p08_candidate_table_path={expected_candidate_table}",
            )
        if not _same_normalized_path(row["scan_manifest_path"], expected_scan_manifest):
            _fail(
                outdir,
                blocks,
                "P07 scan_manifest_path mismatch",
                proteome_shard=key[0],
                target_id=key[1],
                source_path=str(p07_sequence_table),
                notes=f"p07_scan_manifest_path={row['scan_manifest_path']}; p08_scan_manifest_path={expected_scan_manifest}",
            )
        if key in p07_by_key:
            _fail(outdir, blocks, "duplicate P07 join key", proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p07_sequence_id={row['p07_sequence_id']}")
        if row["p07_sequence_id"] in p07_ids:
            _fail(outdir, blocks, "duplicate P07 sequence ID", proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p07_sequence_id={row['p07_sequence_id']}")
        p07_by_key[key] = row
        p07_ids.add(row["p07_sequence_id"])
    statuses: dict[tuple[str, str], dict[str, str]] = {}
    for row in status_rows:
        status_key = (row["tool"], row["fasta_shard"])
        if status_key in statuses:
            _fail(
                outdir,
                blocks,
                "duplicate P07 status key",
                source_path=str(p07_status_table),
                notes=f"tool={status_key[0]}; fasta_shard={status_key[1]}",
            )
        statuses[status_key] = row
    try:
        taxonomy = _load_taxonomy_records(
            [Path(path) for path in taxonomy_paths], taxonomy_source_roles
        )
    except (OSError, ValueError) as error:
        _fail(outdir, blocks, f"taxonomy loading failed: {error}", source_path=", ".join(map(str, taxonomy_paths)))

    registry_by_family = {row["family_category"]: row for row in registry_rows}
    references_by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    seed_reference_by_accession: dict[str, dict[str, str]] = {}
    for record_kind, rows, identifier_field, provenance_table in (("seed", seed_rows, "seed_id", p05_seed_table), ("control", control_rows, "control_id", p05_control_table)):
        for row in rows:
            path = Path(row["sequence_path"])
            source_model = registry_by_family.get(row["family_category"])
            if source_model is None or row["model_sha256"] != source_model["model_sha256"]:
                _fail(outdir, blocks, "P05 reference model SHA-256 mismatch", family_category=row["family_category"], source_path=str(path), notes=f"{identifier_field}={row[identifier_field]}; declared_model_sha256={row['model_sha256']}; registry_model_sha256={(source_model or {}).get('model_sha256', 'missing')}")
            try:
                verified_file_sha256 = _sha256(path)
            except OSError as error:
                _fail(outdir, blocks, f"sequence file unreadable: {path}", family_category=row["family_category"], source_path=str(path), notes=str(error))
            try:
                fasta_records = _read_fasta(path)
            except (OSError, ValueError) as error:
                _fail(outdir, blocks, "reference FASTA malformed", family_category=row["family_category"], source_path=str(path), notes=f"{identifier_field}={row[identifier_field]}; {error}")
            if len(fasta_records) != 1:
                _fail(outdir, blocks, "reference FASTA malformed", family_category=row["family_category"], source_path=str(path), notes=f"{identifier_field}={row[identifier_field]}; expected exactly one FASTA record")
            sequence = next(iter(fasta_records.values()))
            checksum_kind = row.get("source_checksum_kind", "file_sha256")
            verified = _residue_sha256(sequence) if checksum_kind == "residue_sha256" else verified_file_sha256
            if verified != row["sequence_sha256"]:
                _fail(outdir, blocks, "SHA-256 mismatch", family_category=row["family_category"], source_path=str(path), notes=f"{identifier_field}={row[identifier_field]}; checksum_kind={checksum_kind}")
            reference = {
                "family_category": row["family_category"], "record_kind": record_kind, "record_id": row[identifier_field], "source_accession": row.get("source_accession", ""), "control_role": row.get("control_role", ""), "sequence_path": str(path), "sequence_sha256": row["sequence_sha256"], "verified_sha256": verified, "source_checksum_kind": checksum_kind, "model_sha256": row["model_sha256"], "source_model_sha256": row["model_sha256"], "model_provenance": "direct_p05_registry", "model_provenance_source_path": str(provenance_table), "evidence": row.get("evidence", ""), "notes": row.get("notes", ""), "sequence": sequence, "source_accession_or_assembly": row.get("source_accession", "") or row[identifier_field],
            }
            references_by_family[row["family_category"]].append(reference)
            if record_kind == "seed":
                accession = row.get("source_accession", "")
                if accession in seed_reference_by_accession:
                    _fail(outdir, blocks, "duplicate P05 seed source accession", family_category=row["family_category"], source_path=str(path), notes=f"source_accession={accession}")
                seed_reference_by_accession[accession] = reference

    if CORE_FAMILY in selected_families:
        core_model = approved_models[CORE_FAMILY]
        core_seed_path = Path(p05_seed_table).parent / CORE_SEED_REGISTRY_FILENAME
        try:
            core_seed_rows = _read_tsv(core_seed_path, CORE_SEED_FIELDS, "P05 extracellular core seed registry")
        except ValueError as error:
            _fail(outdir, blocks, "missing authoritative core seed registry", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=str(error))
        expected_count = core_model.get("seed_sequence_count", "")
        if (
            len(core_seed_rows) != CORE_SEED_COUNT
            or (expected_count and (not expected_count.isdigit() or int(expected_count) != CORE_SEED_COUNT))
        ):
            _fail(outdir, blocks, "core seed registry count mismatch", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=f"expected_core_seed_count={CORE_SEED_COUNT}; registry_seed_sequence_count={expected_count or 'not_declared'}; core_seed_rows={len(core_seed_rows)}")
        seen_core_accessions: set[str] = set()
        for row in core_seed_rows:
            accession = row["source_accession"]
            if row["family_category"] != CORE_FAMILY or row["model_sha256"] != core_model["model_sha256"]:
                _fail(outdir, blocks, "core seed registry model SHA-256 mismatch", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=f"seed_id={row['seed_id']}; declared_model_sha256={row['model_sha256']}; registry_model_sha256={core_model['model_sha256']}")
            if accession in seen_core_accessions or accession not in seed_reference_by_accession:
                _fail(outdir, blocks, "core seed registry source accession mismatch", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=f"seed_id={row['seed_id']}; source_accession={accession}")
            source = seed_reference_by_accession[accession]
            if Path(row["sequence_path"]) != Path(source["sequence_path"]):
                _fail(outdir, blocks, "core seed registry sequence path mismatch", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=f"seed_id={row['seed_id']}; source_accession={accession}")
            seen_core_accessions.add(accession)
            references_by_family[CORE_FAMILY].append({
                **source, "family_category": CORE_FAMILY, "record_kind": "seed", "record_id": row["seed_id"],
                "model_sha256": core_model["model_sha256"], "source_model_sha256": source["model_sha256"],
                "model_provenance": "derived_core_seed_registry", "model_provenance_source_path": str(core_seed_path),
            })
        cross_family_sources = sorted(
            (
                row
                for family, records in references_by_family.items()
                if family in approved_models and family != CORE_FAMILY
                for row in records
                if row["record_kind"] == "seed"
            ),
            key=lambda row: (row["family_category"], row["source_accession"]),
        )
        if len(cross_family_sources) != CORE_CROSS_FAMILY_CHALLENGE_COUNT:
            _fail(outdir, blocks, "core hard-panel count mismatch", family_category=CORE_FAMILY, source_path=str(p05_seed_table), notes=f"expected_cross_family_challenges={CORE_CROSS_FAMILY_CHALLENGE_COUNT}; observed={len(cross_family_sources)}")
        for source in cross_family_sources:
            references_by_family[CORE_FAMILY].append({
                **source, "family_category": CORE_FAMILY, "record_kind": "control",
                "record_id": f"{CORE_FAMILY}|cross_family_challenge|{source['family_category']}|{source['source_accession']}",
                "control_role": "cross_family_challenge", "model_sha256": core_model["model_sha256"],
                "source_model_sha256": source["model_sha256"], "model_provenance": "derived_core_cross_family_control", "model_provenance_source_path": str(p05_seed_table),
            })
        close_controls_path = Path(p05_seed_table).parent / CORE_CLOSE_CONTROLS_FILENAME
        try:
            close_rows = _read_tsv(close_controls_path, CORE_CLOSE_CONTROL_FIELDS, "P05 extracellular core close controls")
        except ValueError as error:
            _fail(outdir, blocks, "missing authoritative core close controls", family_category=CORE_FAMILY, source_path=str(close_controls_path), notes=str(error))
        if len(close_rows) != CORE_CLOSE_CONTROL_COUNT or len({row["source_accession"] for row in close_rows}) != CORE_CLOSE_CONTROL_COUNT:
            _fail(outdir, blocks, "core hard-panel count mismatch", family_category=CORE_FAMILY, source_path=str(close_controls_path), notes=f"expected_close_controls={CORE_CLOSE_CONTROL_COUNT}; observed={len(close_rows)}")
        for row in close_rows:
            match = re.search(r"residue SHA256 is ([0-9a-f]{64})", row["notes"])
            if match is None:
                _fail(outdir, blocks, "core close-control checksum provenance missing", family_category=CORE_FAMILY, source_path=str(close_controls_path), notes=f"source_accession={row['source_accession']}")
            path = Path(row["sequence_path"])
            try:
                fasta_records = _read_fasta(path)
            except (OSError, ValueError) as error:
                _fail(outdir, blocks, "reference FASTA malformed", family_category=CORE_FAMILY, source_path=str(path), notes=f"source_accession={row['source_accession']}; {error}")
            if len(fasta_records) != 1:
                _fail(outdir, blocks, "reference FASTA malformed", family_category=CORE_FAMILY, source_path=str(path), notes=f"source_accession={row['source_accession']}; expected exactly one FASTA record")
            sequence = next(iter(fasta_records.values()))
            expected_sha256 = match.group(1)
            if _residue_sha256(sequence) != expected_sha256:
                _fail(outdir, blocks, "SHA-256 mismatch", family_category=CORE_FAMILY, source_path=str(path), notes=f"source_accession={row['source_accession']}; checksum_kind=residue_sha256")
            references_by_family[CORE_FAMILY].append({
                "family_category": CORE_FAMILY, "record_kind": "control", "record_id": f"{CORE_FAMILY}|close_non_target_hydrolase|{row['family_category']}|{row['source_accession']}", "source_accession": row["source_accession"], "control_role": "close_non_target_hydrolase", "sequence_path": str(path), "sequence_sha256": expected_sha256, "verified_sha256": expected_sha256, "source_checksum_kind": "residue_sha256", "model_sha256": core_model["model_sha256"], "source_model_sha256": "", "model_provenance": "derived_core_close_control", "model_provenance_source_path": str(close_controls_path), "evidence": row.get("evidence_level", ""), "notes": row["notes"], "sequence": sequence, "source_accession_or_assembly": row["source_accession"],
            })
        core_controls = [row for row in references_by_family[CORE_FAMILY] if row["record_kind"] == "control"]
        if len(core_controls) != CORE_HARD_PANEL_COUNT:
            _fail(outdir, blocks, "core hard-panel count mismatch", family_category=CORE_FAMILY, notes=f"expected_hard_panel={CORE_HARD_PANEL_COUNT}; observed={len(core_controls)}")

    candidate_rows: list[dict[str, str]] = []
    candidate_fasta_records: list[dict[str, str]] = []
    taxonomy_rows: list[dict[str, str]] = []
    candidate_fasta_cache: dict[str, tuple[str, dict[str, str]]] = {}
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
        p07_families = {value.strip() for value in p07["family_categories"].split(";") if value.strip()}
        if family not in p07_families:
            _fail(outdir, blocks, "P07 family_categories missing P06 family", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p07_family_categories={p07['family_categories']}")
        length_text = (p06["target_length"], p07["target_length_from_p06"], p07["sequence_length"])
        if any(not value.isdigit() or int(value) <= 0 for value in length_text):
            _fail(outdir, blocks, "P06/P07 target length is not a positive integer", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p06_target_length={length_text[0]}; p07_target_length_from_p06={length_text[1]}; p07_sequence_length={length_text[2]}")
        expected_length = int(p06["target_length"])
        if len({expected_length, int(p07["target_length_from_p06"]), int(p07["sequence_length"])}) != 1:
            _fail(outdir, blocks, "P06/P07 sequence length mismatch", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table))
        fasta_path = Path(p07["fasta_shard"])
        cached = candidate_fasta_cache.get(str(fasta_path))
        if cached is None:
            try:
                cached = (_sha256(fasta_path), _read_fasta(fasta_path))
            except OSError as error:
                _fail(outdir, blocks, "candidate FASTA unreadable", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(fasta_path), notes=str(error))
            except ValueError as error:
                _fail(outdir, blocks, "candidate FASTA malformed", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(fasta_path), notes=str(error))
            candidate_fasta_cache[str(fasta_path)] = cached
        fasta_sha256, candidate_source_records = cached
        sequence = candidate_source_records.get(p07["p07_sequence_id"])
        if sequence is None:
            _fail(outdir, blocks, "candidate FASTA missing p07_sequence_id", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(fasta_path), notes=f"p07_sequence_id={p07['p07_sequence_id']}")
        if len(sequence) != expected_length:
            _fail(outdir, blocks, "candidate FASTA sequence length mismatch", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(fasta_path), notes=f"p07_sequence_id={p07['p07_sequence_id']}; expected_length={p06['target_length']}; fasta_length={len(sequence)}")
        shard_stem = Path(p07["fasta_shard"]).stem
        required_statuses = {tool: statuses.get((tool, shard_stem)) for tool in REQUIRED_P07_TOOLS}
        status_detail = "; ".join(f"{tool}={(required_statuses[tool] or {}).get('status', 'missing')}" for tool in REQUIRED_P07_TOOLS)
        if any((required_statuses[tool] or {}).get("status") not in {"completed", "skipped_existing"} for tool in REQUIRED_P07_TOOLS):
            _fail(outdir, blocks, "P07 annotation status requirement failed", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["fasta_shard"], notes=f"fasta_shard_stem={shard_stem}; {status_detail}")
        for tool, status_row in required_statuses.items():
            assert status_row is not None
            try:
                input_matches = _same_normalized_path(status_row["input_fasta"], fasta_path)
                _normalized_path(status_row["output_path"])
            except (OSError, ValueError) as error:
                _fail(outdir, blocks, "P07 status path is not parseable", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_status_table), notes=f"tool={tool}; {error}")
            if not input_matches:
                _fail(outdir, blocks, "P07 status input FASTA path mismatch", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=str(p07_status_table), notes=f"tool={tool}; status_input_fasta={status_row['input_fasta']}; p07_fasta_shard={p07['fasta_shard']}")
        assembly_accession = _assembly_accession(p07["source_proteome_path"])
        p03_manifest = p03_manifest_by_accession.get(assembly_accession)
        p03_qc = p03_qc_by_accession.get(assembly_accession)
        if p03_manifest is None or p03_qc is None:
            _fail(outdir, blocks, "P03 provenance absence blocks processing", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["source_proteome_path"], notes=f"assembly_accession={assembly_accession}; manifest={'present' if p03_manifest else 'missing'}; qc={'present' if p03_qc else 'missing'}")
        if Path(p03_manifest["faa_path"]) != Path(p07["source_proteome_path"]) or Path(p03_qc["faa_path"]) != Path(p07["source_proteome_path"]):
            _fail(outdir, blocks, "P03 FAA source path mismatch", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["source_proteome_path"], notes=f"p03_manifest_faa_path={p03_manifest['faa_path']}; p03_qc_faa_path={p03_qc['faa_path']}")
        taxonomy_record = taxonomy.get(assembly_accession)
        if not taxonomy_record:
            _fail(outdir, blocks, "taxonomy absence blocks processing", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["source_proteome_path"])
        lineage = taxonomy_record["lineage"]
        candidate = {
            **{field: p06[field] for field in P06_FIELDS},
            "gtdb_release": expected_gtdb_release, "p07_declared_gtdb_release": p07["gtdb_release"],
            "p06_candidate_table_path": str(expected_candidate_table), "p06_candidate_table_sha256": expected_candidate_table_sha256,
            "p07_declared_candidate_table_path": p07["candidate_table_path"], "p07_declared_candidate_table_sha256": _sha256(Path(p07["candidate_table_path"])),
            "p06_scan_manifest_path": str(expected_scan_manifest), "p06_scan_manifest_sha256": expected_scan_manifest_sha256,
            "p07_declared_scan_manifest_path": p07["scan_manifest_path"], "p07_declared_scan_manifest_sha256": _sha256(Path(p07["scan_manifest_path"])),
            "p03_prediction_manifest_path": str(provenance_paths["p03_prediction_manifest"]), "p03_prediction_manifest_sha256": _sha256(provenance_paths["p03_prediction_manifest"]),
            "p03_prediction_qc_path": str(provenance_paths["p03_prediction_qc"]), "p03_prediction_qc_sha256": _sha256(provenance_paths["p03_prediction_qc"]), "p03_faa_path": p03_manifest["faa_path"],
            "p07_sequence_id": p07["p07_sequence_id"], "source_proteome_path": p07["source_proteome_path"], "fasta_shard": p07["fasta_shard"],
            "p07_annotation_status": ";".join(sorted({required_statuses[tool]["status"] for tool in REQUIRED_P07_TOOLS})),
            "p07_annotation_status_by_tool": ";".join(f"{tool}={required_statuses[tool]['status']}" for tool in REQUIRED_P07_TOOLS),
            "p07_annotation_status_table_path": str(p07_status_table), "p07_annotation_status_table_sha256": _sha256(Path(p07_status_table)),
            "p07_annotation_output_paths": ";".join(f"{tool}={required_statuses[tool]['output_path']}" for tool in REQUIRED_P07_TOOLS),
            "p07_annotation_input_fastas": ";".join(f"{tool}={required_statuses[tool]['input_fasta']}" for tool in REQUIRED_P07_TOOLS),
            "assembly_accession": assembly_accession, "taxonomy_lineage": lineage, "model_sha256": approved_models[family]["model_sha256"], "evidence_boundary": EVIDENCE_BOUNDARY,
        }
        candidate_rows.append(candidate)
        candidate_fasta_records.append({
            "family_category": family, "record_kind": "candidate", "record_id": p07["p07_sequence_id"],
            "sequence": sequence, "source_path": str(fasta_path), "source_sha256": fasta_sha256,
            "model_sha256": approved_models[family]["model_sha256"], "p06_proteome_shard": p06["proteome_shard"],
            "p06_target_id": p06["target_id"], "p07_sequence_id": p07["p07_sequence_id"],
            "source_accession_or_assembly": assembly_accession, "evidence_or_control_role": p06["tier"],
            "is_gtdb_candidate": "yes", "source_checksum_kind": "file_sha256",
            "source_model_sha256": approved_models[family]["model_sha256"], "model_provenance": "candidate_from_approved_model", "model_provenance_source_path": str(p05_model_registry),
        })
        taxonomy_rows.append({
            "family_category": family,
            "proteome_shard": p06["proteome_shard"],
            "target_id": p06["target_id"],
            "assembly_accession": assembly_accession,
            "taxonomy_lineage": lineage,
            "taxonomy_source_role": taxonomy_record["source_role"],
            "taxonomy_source_path": taxonomy_record["source_path"],
            "taxonomy_source_sha256": taxonomy_record["source_sha256"],
        })

    sort_key = lambda row: (row["family_category"], row["proteome_shard"], row["target_id"])
    candidate_rows.sort(key=sort_key)
    taxonomy_rows.sort(key=sort_key)
    reference_rows = sorted((row for family in {row["family_category"] for row in selected} for row in references_by_family[family]), key=lambda row: (row["family_category"], row["record_kind"], row["record_id"]))
    outputs = {
        "candidate_manifest": outdir / "manifests" / "p08_candidate_manifest.tsv",
        "taxonomy_join": outdir / "gtdb_mapping" / "p08_taxonomy_join.tsv",
        "family_reference_manifest": outdir / "manifests" / "p08_family_reference_manifest.tsv",
        "preparation_summary": outdir / "manifests" / "p08_preparation_summary.tsv",
        "family_input_manifest": outdir / "manifests" / "p08_family_input_manifest.tsv",
        "phylogeny_command_manifest": outdir / "manifests" / "p08_phylogeny_command_manifest.tsv",
        "input_provenance": outdir / "manifests" / "p08_input_provenance.tsv",
    }
    family_input_rows: list[dict[str, str]] = []
    command_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    kind_order = {"candidate": 0, "seed": 1, "control": 2}
    for family in sorted({row["family_category"] for row in candidate_rows}):
        candidate_inputs = [row for row in candidate_fasta_records if row["family_category"] == family]
        reference_inputs = []
        for reference in references_by_family[family]:
            reference_inputs.append({
                "family_category": family, "record_kind": reference["record_kind"], "record_id": reference["record_id"],
                "sequence": reference["sequence"], "source_path": reference["sequence_path"],
                "source_sha256": reference["verified_sha256"], "source_checksum_kind": reference["source_checksum_kind"],
                "model_sha256": reference["model_sha256"], "source_model_sha256": reference["source_model_sha256"],
                "model_provenance": reference["model_provenance"], "model_provenance_source_path": reference["model_provenance_source_path"],
                "p06_proteome_shard": "", "p06_target_id": "", "p07_sequence_id": "",
                "source_accession_or_assembly": reference["source_accession_or_assembly"],
                "evidence_or_control_role": reference["control_role"] or reference["evidence"], "is_gtdb_candidate": "no",
            })
        family_inputs = sorted(candidate_inputs + reference_inputs, key=lambda row: (kind_order[row["record_kind"]], row["record_id"], row["source_path"]))
        fasta_path = outdir / "family_fastas" / f"{family}.faa"
        _write_fasta(fasta_path, family_inputs)
        input_sha256 = _sha256(fasta_path)
        for input_row in family_inputs:
            record_identity = "|".join((input_row["record_kind"], family, input_row["p06_proteome_shard"], input_row["p06_target_id"], input_row["record_id"]))
            family_input_rows.append({
                "family_category": family, "record_kind": input_row["record_kind"], "record_identity": record_identity,
                "input_fasta_path": str(fasta_path), "input_sha256": input_sha256, "source_path": input_row["source_path"],
                "source_sha256": input_row["source_sha256"], "source_checksum_kind": input_row["source_checksum_kind"],
                "model_sha256": input_row["model_sha256"], "source_model_sha256": input_row["source_model_sha256"],
                "model_provenance": input_row["model_provenance"], "model_provenance_source_path": input_row["model_provenance_source_path"],
                "p06_proteome_shard": input_row["p06_proteome_shard"], "p06_target_id": input_row["p06_target_id"],
                "p07_sequence_id": input_row["p07_sequence_id"], "source_accession_or_assembly": input_row["source_accession_or_assembly"],
                "evidence_or_control_role": input_row["evidence_or_control_role"], "is_gtdb_candidate": input_row["is_gtdb_candidate"],
                "evidence_boundary": EVIDENCE_BOUNDARY,
            })
        candidate_count = len(candidate_inputs)
        templates = planned_command_templates(
            candidate_count,
            mafft_exe=mafft_exe,
            iqtree_exe=iqtree_exe,
            fasttree_exe=fasttree_exe,
        )
        alignment_fasta = outdir / "alignments" / f"{family}.aligned.faa"
        representative_input = outdir / "review" / f"{family}.representative_input.faa"
        representative_mapping = outdir / "review" / f"{family}.representative_selection_mapping.tsv"
        fasttree_tree = outdir / "trees" / f"{family}.fasttree.nwk"
        iqtree_prefix = outdir / "trees" / family
        command_rows.append({
            "family_category": family, "command_status": "planned_not_run", "input_fasta_path": str(fasta_path),
            "input_sha256": input_sha256, "candidate_input_record_count": str(candidate_count),
            "total_input_record_count": str(len(family_inputs)), "route": templates["route"],
            "alignment_fasta_path": str(alignment_fasta), "representative_input_fasta_path": str(representative_input),
            "fasttree_tree_path": str(fasttree_tree), "iqtree_prefix": str(iqtree_prefix),
            "representative_plan": templates["representative_plan"], "mafft_template": templates["mafft_template"],
            "representative_selection_algorithm": templates["representative_selection_algorithm"],
            "representative_selection_algorithm_version": templates["representative_selection_algorithm_version"],
            "representative_selection_parameters": templates["representative_selection_parameters"],
            "representative_selection_mapping_path": str(representative_mapping) if templates["representative_selection_algorithm"] else "",
            "representative_selection_mapping_sha256": templates["representative_selection_mapping_sha256"],
            "representative_selection_mapping_record_count": templates["representative_selection_mapping_record_count"],
            "representative_materialization_status": templates["representative_materialization_status"],
            "fasttree_template": templates["fasttree_template"], "iqtree2_template": templates["iqtree2_template"],
            "iqtree2_annotation": "requires_independent_subset_and_outgroup_approval", "rooting_policy": ROOTING_POLICY,
            "evidence_boundary": EVIDENCE_BOUNDARY,
        })
        summary_rows.append({
            "family_category": family, "candidate_count": str(candidate_count), "route": templates["route"],
            "family_fasta_path": str(fasta_path), "family_fasta_sha256": input_sha256,
            "total_fasta_record_count": str(len(family_inputs)),
        })
    _write_tsv(outputs["candidate_manifest"], tuple(candidate_rows[0]), candidate_rows)
    _write_tsv(outputs["taxonomy_join"], TAXONOMY_JOIN_FIELDS, taxonomy_rows)
    _write_tsv(outputs["family_reference_manifest"], ("family_category", "record_kind", "record_id", "source_accession", "control_role", "sequence_path", "sequence_sha256", "verified_sha256", "evidence", "notes"), reference_rows)
    _write_tsv(outputs["preparation_summary"], SUMMARY_FIELDS, summary_rows)
    _write_tsv(outputs["family_input_manifest"], FAMILY_INPUT_FIELDS, sorted(family_input_rows, key=lambda row: (row["family_category"], kind_order[row["record_kind"]], row["record_identity"])))
    _write_tsv(outputs["phylogeny_command_manifest"], COMMAND_FIELDS, command_rows)
    source_roles = tuple(taxonomy_source_roles or ())
    primary_inputs = (
        ("candidate_table", Path(p06_candidate_table), "P06_candidate_input"),
        ("p06_scan_manifest", provenance_paths["p06_scan_manifest"], "P06_accepted_scan_contract"),
        ("p07_sequence_manifest", Path(p07_sequence_table), "P07_sequence_input"),
        ("p07_status_table", Path(p07_status_table), "P07_annotation_status_input"),
        ("p03_prediction_manifest", provenance_paths["p03_prediction_manifest"], "P03_prediction_and_FAA_source_contract"),
        ("p03_prediction_qc", provenance_paths["p03_prediction_qc"], "P03_prediction_QC_contract"),
        ("model_registry", Path(p05_model_registry), "P05_approved_model_gate"),
        ("seed_registry", Path(p05_seed_table), "P05_accessioned_seed_provenance"),
        ("control_panel", Path(p05_control_table), "P05_accessioned_control_provenance"),
    )
    taxonomy_inputs = tuple(
        (source_roles[index] if source_roles else f"taxonomy_input_{index + 1}", Path(path), "GTDB_taxonomy_mapping")
        for index, path in enumerate(taxonomy_paths)
    )
    tree_inputs = tuple(
        (role, Path(path), "provenance_preflight_only_no_topology_read")
        for role, path in (additional_provenance_inputs or {}).items()
    )
    provenance_rows = [
        {"input_role": role, "input_path": str(path), "input_sha256": _sha256(path), "input_usage": usage}
        for role, path, usage in (*primary_inputs, *taxonomy_inputs, *tree_inputs)
    ]
    provenance_rows.append({
        "input_role": "gtdb_release",
        "input_path": gtdb_release.strip(),
        "input_sha256": "not_a_file_declaration",
        "input_usage": "explicit_GTDB_release_declaration_not_inferred",
    })
    _write_tsv(outputs["input_provenance"], INPUT_PROVENANCE_FIELDS, sorted(provenance_rows, key=lambda row: row["input_role"]))
    return outputs


def _require_existing_provenance_path(path: Path, label: str) -> Path:
    """Require a local provenance/preflight input without opening or inferring from it."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{label} must be an existing nonempty provenance/preflight input: {path}")
    return path


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the planning-only P08 preparation CLI parser."""
    parser = argparse.ArgumentParser(
        description="Prepare P08 phylogeny manifests without running alignment or tree inference."
    )
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--gtdb-release", required=True, help="Explicit GTDB release declaration; it is recorded and never inferred.")
    parser.add_argument("--p06-scan-manifest", type=Path, required=True)
    parser.add_argument("--p03-prediction-manifest", type=Path, required=True)
    parser.add_argument("--p03-prediction-qc", type=Path, required=True)
    parser.add_argument("--p07-sequence-manifest", type=Path, required=True)
    parser.add_argument("--p07-status-table", type=Path, required=True)
    parser.add_argument("--model-registry", type=Path, required=True)
    parser.add_argument("--seed-registry", type=Path, required=True)
    parser.add_argument("--control-panel", type=Path, required=True)
    parser.add_argument("--bac120-taxonomy", type=Path, required=True)
    parser.add_argument("--ar53-taxonomy", type=Path, required=True)
    parser.add_argument("--bac120-tree", type=Path, required=True)
    parser.add_argument("--ar53-tree", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--include-tier", choices=("High-confidence", "Review"), action="append")
    parser.add_argument("--mafft-exe", default="mafft")
    parser.add_argument("--iqtree-exe", default="iqtree2")
    parser.add_argument("--fasttree-exe", default="FastTree")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Write immutable P08 planning manifests; external tools are never invoked."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        bac120_tree = _require_existing_provenance_path(args.bac120_tree, "--bac120-tree")
        ar53_tree = _require_existing_provenance_path(args.ar53_tree, "--ar53-tree")
        outputs = prepare_p08_inputs(
            p06_candidate_table=args.candidate_table,
            gtdb_release=args.gtdb_release,
            p06_scan_manifest=args.p06_scan_manifest,
            p03_prediction_manifest=args.p03_prediction_manifest,
            p03_prediction_qc=args.p03_prediction_qc,
            p07_sequence_table=args.p07_sequence_manifest,
            p07_status_table=args.p07_status_table,
            p05_model_registry=args.model_registry,
            p05_seed_table=args.seed_registry,
            p05_control_table=args.control_panel,
            taxonomy_paths=(args.bac120_taxonomy, args.ar53_taxonomy),
            taxonomy_source_roles=("bac120_taxonomy", "ar53_taxonomy"),
            additional_provenance_inputs={"bac120_tree": bac120_tree, "ar53_tree": ar53_tree},
            outdir=args.outdir,
            include_tiers=tuple(args.include_tier or DEFAULT_INCLUDE_TIERS),
            mafft_exe=args.mafft_exe,
            iqtree_exe=args.iqtree_exe,
            fasttree_exe=args.fasttree_exe,
        )
    except ValueError as error:
        parser.error(str(error))
    print("P08 preparation complete: planned_not_run only; no MAFFT, IQ-TREE, or FastTree command was executed.")
    print(f"Bac120 tree provenance/preflight input: {bac120_tree}")
    print(f"Ar53 tree provenance/preflight input: {ar53_tree}")
    for name, path in sorted(outputs.items()):
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
