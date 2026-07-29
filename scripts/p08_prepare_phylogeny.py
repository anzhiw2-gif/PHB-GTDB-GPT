"""Prepare auditable P08 phylogeny inputs without running phylogeny tools.

This module joins accepted P06 candidates to P07 sequence/annotation records,
GTDB taxonomy, and checksum-verified P05 reference records.  It deliberately
does not execute alignment, tree, or other external phylogeny tools.
"""

from __future__ import annotations

import csv
import hashlib
import re
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
SEED_FIELDS = ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path", "sequence_sha256")
CONTROL_FIELDS = ("family_category", "model_sha256", "control_id", "control_role", "sequence_path", "sequence_sha256")
CORE_FAMILY = "extracellular_pha_depolymerase_core"
CORE_SEED_REGISTRY_FILENAME = "p05_extracellular_core_seed_registry.tsv"
CORE_CLOSE_CONTROLS_FILENAME = "p05_extracellular_core_close_controls.tsv"
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
    "mafft_template", "fasttree_template", "iqtree2_template", "iqtree2_annotation",
    "rooting_policy", "evidence_boundary",
)
SUMMARY_FIELDS = (
    "family_category", "candidate_count", "route", "family_fasta_path", "family_fasta_sha256",
    "total_fasta_record_count",
)
ROOTING_POLICY = "explicit_accessioned_outgroup_required; otherwise midpoint_display_only"
IQTREE2_TEMPLATE = "iqtree2 -s {alignment_fasta} -m TEST -B 1000 --prefix {iqtree_prefix}"


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


def planned_command_templates(candidate_count: int) -> dict[str, str]:
    """Return unexecuted command templates appropriate to a candidate count."""
    route = route_family_size(candidate_count)
    templates = {
        "route": route,
        "mafft_template": "",
        "fasttree_template": "",
        "representative_plan": "",
        "iqtree2_template": IQTREE2_TEMPLATE,
    }
    if route == "mafft_linsi_then_review":
        templates["mafft_template"] = "mafft --localpair --maxiterate 1000 --thread {threads} --inputorder {input_fasta} > {alignment_fasta}"
    elif route == "mafft_auto_then_review":
        templates["mafft_template"] = "mafft --auto --thread {threads} --inputorder {input_fasta} > {alignment_fasta}"
    else:
        templates["representative_plan"] = "deterministic representative subset required before exploratory FastTree; no representative selection executed"
        templates["fasttree_template"] = "FastTree -lg {representative_alignment_fasta} > {fasttree_tree}"
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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

    p07_by_key: dict[tuple[str, str], dict[str, str]] = {}
    p07_ids: set[str] = set()
    for row in p07_rows:
        key = (row["proteome_shard"], row["target_id"])
        if key in p07_by_key:
            _fail(outdir, blocks, "duplicate P07 join key", proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p07_sequence_id={row['p07_sequence_id']}")
        if row["p07_sequence_id"] in p07_ids:
            _fail(outdir, blocks, "duplicate P07 sequence ID", proteome_shard=key[0], target_id=key[1], source_path=str(p07_sequence_table), notes=f"p07_sequence_id={row['p07_sequence_id']}")
        p07_by_key[key] = row
        p07_ids.add(row["p07_sequence_id"])
    statuses = {(row["tool"], Path(row["fasta_shard"]).stem): row["status"] for row in status_rows}
    try:
        taxonomy = load_taxonomy([Path(path) for path in taxonomy_paths])
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
            checksum_kind = "residue_sha256" if record_kind == "control" and "hard_negative" in row else "file_sha256"
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

    selected_families = {row["family_category"] for row in selected}
    if CORE_FAMILY in selected_families and not any(row["record_kind"] == "seed" for row in references_by_family[CORE_FAMILY]):
        core_model = approved_models[CORE_FAMILY]
        core_seed_path = Path(p05_seed_table).parent / CORE_SEED_REGISTRY_FILENAME
        try:
            core_seed_rows = _read_tsv(core_seed_path, CORE_SEED_FIELDS, "P05 extracellular core seed registry")
        except ValueError as error:
            _fail(outdir, blocks, "missing authoritative core seed registry", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=str(error))
        expected_count = core_model.get("seed_sequence_count", "")
        if expected_count and (not expected_count.isdigit() or len(core_seed_rows) != int(expected_count)):
            _fail(outdir, blocks, "core seed registry count mismatch", family_category=CORE_FAMILY, source_path=str(core_seed_path), notes=f"registry_seed_sequence_count={expected_count}; core_seed_rows={len(core_seed_rows)}")
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
        if not any(row["record_kind"] == "control" for row in references_by_family[CORE_FAMILY]):
            for source in sorted((row for family, records in references_by_family.items() if family in approved_models and family != CORE_FAMILY for row in records if row["record_kind"] == "seed"), key=lambda row: (row["family_category"], row["source_accession"])):
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
        status_detail = "; ".join(f"{tool}={statuses.get((tool, shard_stem), 'missing')}" for tool in REQUIRED_P07_TOOLS)
        if any(statuses.get((tool, shard_stem)) not in {"completed", "skipped_existing"} for tool in REQUIRED_P07_TOOLS):
            _fail(outdir, blocks, "P07 annotation status requirement failed", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["fasta_shard"], notes=f"fasta_shard_stem={shard_stem}; {status_detail}")
        assembly_accession = _assembly_accession(p07["source_proteome_path"])
        lineage = taxonomy.get(assembly_accession)
        if not lineage:
            _fail(outdir, blocks, "taxonomy absence blocks processing", family_category=family, proteome_shard=key[0], target_id=key[1], source_path=p07["source_proteome_path"])
        candidate = {
            **{field: p06[field] for field in P06_FIELDS},
            "p07_sequence_id": p07["p07_sequence_id"], "source_proteome_path": p07["source_proteome_path"], "fasta_shard": p07["fasta_shard"], "p07_annotation_status": "completed", "assembly_accession": assembly_accession, "taxonomy_lineage": lineage, "model_sha256": approved_models[family]["model_sha256"], "evidence_boundary": EVIDENCE_BOUNDARY,
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
        taxonomy_rows.append({"family_category": family, "proteome_shard": p06["proteome_shard"], "target_id": p06["target_id"], "assembly_accession": assembly_accession, "taxonomy_lineage": lineage})

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
        templates = planned_command_templates(candidate_count)
        alignment_fasta = outdir / "alignments" / f"{family}.aligned.faa"
        representative_input = outdir / "review" / f"{family}.representative_input.faa"
        fasttree_tree = outdir / "trees" / f"{family}.fasttree.nwk"
        iqtree_prefix = outdir / "trees" / family
        command_rows.append({
            "family_category": family, "command_status": "planned_not_run", "input_fasta_path": str(fasta_path),
            "input_sha256": input_sha256, "candidate_input_record_count": str(candidate_count),
            "total_input_record_count": str(len(family_inputs)), "route": templates["route"],
            "alignment_fasta_path": str(alignment_fasta), "representative_input_fasta_path": str(representative_input),
            "fasttree_tree_path": str(fasttree_tree), "iqtree_prefix": str(iqtree_prefix),
            "representative_plan": templates["representative_plan"], "mafft_template": templates["mafft_template"],
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
    _write_tsv(outputs["taxonomy_join"], ("family_category", "proteome_shard", "target_id", "assembly_accession", "taxonomy_lineage"), taxonomy_rows)
    _write_tsv(outputs["family_reference_manifest"], ("family_category", "record_kind", "record_id", "source_accession", "control_role", "sequence_path", "sequence_sha256", "verified_sha256", "evidence", "notes"), reference_rows)
    _write_tsv(outputs["preparation_summary"], SUMMARY_FIELDS, summary_rows)
    _write_tsv(outputs["family_input_manifest"], FAMILY_INPUT_FIELDS, sorted(family_input_rows, key=lambda row: (row["family_category"], kind_order[row["record_kind"]], row["record_identity"])))
    _write_tsv(outputs["phylogeny_command_manifest"], COMMAND_FIELDS, command_rows)
    return outputs
