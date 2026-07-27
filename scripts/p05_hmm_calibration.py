"""Materialize accessioned P05 HMM calibration controls.

The panel challenges every model with seeds from the other active families and
keeps target-family boundary records as observations. Boundary records are not
treated as verified negatives, so they never determine a rejection threshold.
All outputs are sequence-evidence records and do not establish a degradation
phenotype for any homologous sequence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


CONTROL_PANEL_FILENAME = "p05_hmm_calibration_control_panel.tsv"
CALIBRATION_COMMAND_MANIFEST_FILENAME = "p05_hmm_calibration_command_manifest.tsv"
CALIBRATION_COMMAND_SUMMARY_FILENAME = "p05_hmm_calibration_command_summary.tsv"
LEAVE_ONE_OUT_COMMAND_MANIFEST_FILENAME = "p05_hmm_leave_one_out_command_manifest.tsv"
LEAVE_ONE_OUT_COMMAND_SUMMARY_FILENAME = "p05_hmm_leave_one_out_command_summary.tsv"
LEAVE_ONE_OUT_RESULT_FILENAME = "p05_hmm_leave_one_out_positive_results.tsv"
CONTROL_SMOKE_RESULT_FILENAME = "p05_hmm_control_smoke_results.tsv"
CALIBRATION_DECISION_FILENAME = "p05_hmm_calibration_decision_summary.tsv"
REFERENCE_REQUIRED_FIELDS = (
    "family_category",
    "source_accession",
    "organism",
    "taxonomic_domain",
    "evidence_level",
    "profile_seed_status",
    "sequence_path",
    "source_database",
    "source_release",
    "source_version",
    "retrieval_date",
    "source_url",
    "doi",
    "pmid",
    "pmcid",
    "literature_support_scope",
    "notes",
)
MODEL_REGISTRY_REQUIRED_FIELDS = ("family_category", "model_sha256")
SEED_REGISTRY_REQUIRED_FIELDS = ("family_category", "source_accession", "model_sha256")
CONTROL_PANEL_FIELDNAMES = (
    "family_category",
    "model_sha256",
    "control_id",
    "control_role",
    "hard_negative",
    "expected_outcome",
    "control_family_category",
    "source_accession",
    "organism",
    "taxonomic_domain",
    "evidence_level",
    "profile_seed_status",
    "sequence_path",
    "sequence_sha256",
    "source_database",
    "source_release",
    "source_version",
    "retrieval_date",
    "source_url",
    "doi",
    "pmid",
    "pmcid",
    "literature_support_scope",
    "architecture_rationale",
    "notes",
)
CALIBRATION_COMMAND_FIELDNAMES = (
    "family_category",
    "model_sha256",
    "hmm_path",
    "target_fasta_path",
    "target_fasta_sha256",
    "domtblout_path",
    "main_output_path",
    "control_count",
    "hard_challenge_count",
    "boundary_observation_count",
    "command",
    "command_status",
    "notes",
)
SUMMARY_FIELDNAMES = ("kind", "name", "count")
COMMAND_STATUS = "planned_not_run"
MINIMUM_LEAVE_ONE_OUT_SEEDS = 4


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sequence_residue_sha256(path: Path) -> str:
    """Hash normalized amino-acid residues, independent of FASTA formatting."""

    return hashlib.sha256(_read_single_fasta_sequence(path).encode("ascii")).hexdigest()


def _load_tsv(path: Path, required_fields: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required calibration input is missing: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing = [field for field in required_fields if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def _load_model_hashes(path: Path) -> dict[str, str]:
    model_hashes: dict[str, str] = {}
    for line_number, row in enumerate(_load_tsv(path, MODEL_REGISTRY_REQUIRED_FIELDS), start=2):
        family = row["family_category"]
        model_hash = row["model_sha256"].lower()
        if not family or not model_hash:
            raise ValueError(f"{path}:{line_number} has an empty family_category or model_sha256")
        if family in model_hashes:
            raise ValueError(f"{path}:{line_number} has duplicate model family {family!r}")
        model_hashes[family] = model_hash
    if len(model_hashes) < 2:
        raise ValueError("Calibration requires at least two active model families for cross-family challenges")
    return model_hashes


def _load_profile_seeds(path: Path, active_families: set[str]) -> dict[str, list[str]]:
    seeds: dict[str, list[str]] = {family: [] for family in active_families}
    for line_number, row in enumerate(_load_tsv(path, SEED_REGISTRY_REQUIRED_FIELDS), start=2):
        family = row["family_category"]
        accession = row["source_accession"]
        if family not in active_families:
            continue
        if not accession:
            raise ValueError(f"{path}:{line_number} has an empty source_accession")
        if accession in seeds[family]:
            raise ValueError(f"{path}:{line_number} has duplicate seed accession {accession!r} for {family}")
        seeds[family].append(accession)
    missing = sorted(family for family, accessions in seeds.items() if not accessions)
    if missing:
        raise ValueError(f"Seed registry has no current profile seeds for: {', '.join(missing)}")
    return {family: sorted(accessions) for family, accessions in seeds.items()}


def _reference_index(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    references: dict[tuple[str, str], dict[str, str]] = {}
    for line_number, row in enumerate(_load_tsv(path, REFERENCE_REQUIRED_FIELDS), start=2):
        key = (row["family_category"], row["source_accession"])
        if not all(key):
            raise ValueError(f"{path}:{line_number} has an empty family_category or source_accession")
        if key in references:
            raise ValueError(f"{path}:{line_number} duplicates reference accession {key!r}")
        references[key] = row
    return references


def _control_row(
    *,
    family: str,
    model_sha256: str,
    reference: dict[str, str],
    control_role: str,
    hard_negative: str,
    expected_outcome: str,
    architecture_rationale: str,
) -> dict[str, str]:
    sequence_path = Path(reference["sequence_path"])
    if not sequence_path.is_file():
        raise FileNotFoundError(
            f"Calibration control {reference['source_accession']!r} for {family} is missing its reference FASTA: {sequence_path}"
        )
    accession = reference["source_accession"]
    control_family = reference["family_category"]
    return {
        "family_category": family,
        "model_sha256": model_sha256,
        "control_id": f"{family}|{control_role}|{control_family}|{accession}",
        "control_role": control_role,
        "hard_negative": hard_negative,
        "expected_outcome": expected_outcome,
        "control_family_category": control_family,
        "source_accession": accession,
        "organism": reference["organism"],
        "taxonomic_domain": reference["taxonomic_domain"],
        "evidence_level": reference["evidence_level"],
        "profile_seed_status": reference["profile_seed_status"],
        "sequence_path": sequence_path.as_posix(),
        "sequence_sha256": _sequence_residue_sha256(sequence_path),
        "source_database": reference["source_database"],
        "source_release": reference["source_release"],
        "source_version": reference["source_version"],
        "retrieval_date": reference["retrieval_date"],
        "source_url": reference["source_url"],
        "doi": reference["doi"],
        "pmid": reference["pmid"],
        "pmcid": reference["pmcid"],
        "literature_support_scope": reference["literature_support_scope"],
        "architecture_rationale": architecture_rationale,
        "notes": reference["notes"],
    }


def build_control_panel(
    reference_manifest: Path,
    seed_registry: Path,
    model_registry: Path,
) -> list[dict[str, str]]:
    """Build deterministic hard challenges and report-only boundary observations.

    A hard challenge is a seed from another active family, where the model must
    not meet the final calibrated threshold. A target-family boundary record is
    biologically ambiguous by construction, so it is reported but cannot set a
    rejection threshold.
    """

    model_hashes = _load_model_hashes(model_registry)
    active_families = set(model_hashes)
    profile_seeds = _load_profile_seeds(seed_registry, active_families)
    references = _reference_index(reference_manifest)
    rows: list[dict[str, str]] = []

    for family in sorted(active_families):
        for control_family in sorted(active_families - {family}):
            for accession in profile_seeds[control_family]:
                reference = references.get((control_family, accession))
                if reference is None:
                    raise ValueError(
                        f"Seed registry accession {accession!r} for {control_family} is absent from {reference_manifest}"
                    )
                rows.append(
                    _control_row(
                        family=family,
                        model_sha256=model_hashes[family],
                        reference=reference,
                        control_role="cross_family_challenge",
                        hard_negative="yes",
                        expected_outcome="must_fail_threshold",
                        architecture_rationale=(
                            f"Current profile seed for distinct active family {control_family}; used to test "
                            "cross-family HMM separation, not as a phenotype-negative assertion."
                        ),
                    )
                )

        for (reference_family, _), reference in sorted(references.items()):
            if reference_family != family or reference["profile_seed_status"] != "boundary_candidate":
                continue
            rows.append(
                _control_row(
                    family=family,
                    model_sha256=model_hashes[family],
                    reference=reference,
                    control_role="boundary_observation",
                    hard_negative="no",
                    expected_outcome="report_only",
                    architecture_rationale=(
                        "Excluded from this profile by the accession-level seed decision; retained to expose "
                        "architecture or evidence ambiguity and never used to set a rejection threshold."
                    ),
                )
            )

    hard_challenge_families = {row["family_category"] for row in rows if row["hard_negative"] == "yes"}
    missing_challenges = sorted(active_families - hard_challenge_families)
    if missing_challenges:
        raise ValueError(f"Calibration panel lacks a cross-family challenge for: {', '.join(missing_challenges)}")
    return sorted(
        rows,
        key=lambda row: (row["family_category"], row["control_role"], row["control_family_category"], row["source_accession"]),
    )


def write_control_panel(path: Path, rows: list[dict[str, str]]) -> Path:
    """Write the compact, Git-trackable calibration panel deterministically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONTROL_PANEL_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in CONTROL_PANEL_FIELDNAMES} for row in rows)
    return path


def _load_models_with_paths(path: Path) -> dict[str, dict[str, str]]:
    rows = _load_tsv(path, MODEL_REGISTRY_REQUIRED_FIELDS + ("model_path",))
    models: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        family = row["family_category"]
        model_path = Path(row["model_path"])
        expected_hash = row["model_sha256"].lower()
        if family in models:
            raise ValueError(f"{path}:{line_number} has duplicate model family {family!r}")
        if not model_path.is_file():
            raise FileNotFoundError(f"{path}:{line_number} model is missing: {model_path}")
        observed_hash = _sha256(model_path)
        if observed_hash != expected_hash:
            raise ValueError(
                f"{path}:{line_number} model SHA256 mismatch for {model_path}: expected {expected_hash}, observed {observed_hash}"
            )
        normalized = dict(row)
        normalized["model_path"] = model_path.as_posix()
        normalized["model_sha256"] = expected_hash
        models[family] = normalized
    return models


def _load_control_panel(path: Path) -> list[dict[str, str]]:
    required = (
        "family_category",
        "model_sha256",
        "control_id",
        "control_role",
        "hard_negative",
        "expected_outcome",
        "source_accession",
        "sequence_path",
        "sequence_sha256",
    )
    rows = _load_tsv(path, required)
    seen_control_ids: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        key = (row["family_category"], row["control_id"])
        if key in seen_control_ids:
            raise ValueError(f"{path}:{line_number} has duplicate calibration control {key!r}")
        seen_control_ids.add(key)
        if row["hard_negative"] == "yes" and row["expected_outcome"] != "must_fail_threshold":
            raise ValueError(f"{path}:{line_number} hard challenge must have expected_outcome=must_fail_threshold")
        if row["hard_negative"] == "no" and row["expected_outcome"] != "report_only":
            raise ValueError(f"{path}:{line_number} boundary observation must have expected_outcome=report_only")
        sequence_path = Path(row["sequence_path"])
        if not sequence_path.is_file():
            raise FileNotFoundError(f"{path}:{line_number} control sequence is missing: {sequence_path}")
        observed_hash = _sequence_residue_sha256(sequence_path)
        if observed_hash != row["sequence_sha256"].lower():
            raise ValueError(
                f"{path}:{line_number} control sequence SHA256 mismatch for {sequence_path}: "
                f"expected {row['sequence_sha256']}, observed {observed_hash}"
            )
    return rows


def _read_single_fasta_sequence(path: Path) -> str:
    headers = 0
    sequence_lines: list[str] = []
    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            headers += 1
            continue
        sequence_lines.append(line)
    sequence = "".join(sequence_lines)
    if headers != 1 or not sequence:
        raise ValueError(f"{path} must contain exactly one nonempty FASTA record for calibration")
    return sequence


def _write_control_fasta(path: Path, rows: list[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            sequence = _read_single_fasta_sequence(Path(row["sequence_path"]))
            handle.write(f">{row['control_id']}\n")
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")
    return path


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _hmmsearch_command(hmm_path: Path, target_fasta_path: Path, domtblout_path: Path, main_output_path: Path) -> str:
    return (
        "hmmsearch --noali --acc --seed 42 --cpu 1 "
        f"--domtblout {_shell_quote(domtblout_path.as_posix())} "
        f"{_shell_quote(hmm_path.as_posix())} {_shell_quote(target_fasta_path.as_posix())} "
        f"> {_shell_quote(main_output_path.as_posix())}"
    )


def _write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fieldnames} for row in rows)
    return path


def build_calibration_command_manifest(
    control_panel: Path,
    model_registry: Path,
    calibration_dir: Path,
) -> dict[str, Path]:
    """Materialize checksum-locked controls and deterministic non-executed commands."""

    models = _load_models_with_paths(model_registry)
    panel_rows = _load_control_panel(control_panel)
    rows_by_family: dict[str, list[dict[str, str]]] = {family: [] for family in models}
    for panel_row in panel_rows:
        family = panel_row["family_category"]
        model = models.get(family)
        if model is None:
            raise ValueError(f"{control_panel} has a control for model absent from registry: {family}")
        if panel_row["model_sha256"].lower() != model["model_sha256"]:
            raise ValueError(f"{control_panel} has a model SHA256 mismatch for {family}")
        rows_by_family[family].append(panel_row)

    target_dir = calibration_dir / "targets"
    raw_dir = calibration_dir / "raw_domtblout"
    log_dir = calibration_dir / "hmmer_logs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    command_rows: list[dict[str, str]] = []
    for family, model in sorted(models.items()):
        family_rows = sorted(rows_by_family[family], key=lambda row: row["control_id"])
        if not family_rows:
            raise ValueError(f"{control_panel} has no controls for {family}")
        hard_count = sum(row["hard_negative"] == "yes" for row in family_rows)
        if hard_count == 0:
            raise ValueError(f"{control_panel} has no cross-family hard challenge for {family}")
        target_fasta_path = _write_control_fasta(target_dir / f"{family}.controls.faa", family_rows)
        domtblout_path = raw_dir / f"{family}.controls.domtblout"
        main_output_path = log_dir / f"{family}.controls.txt"
        hmm_path = Path(model["model_path"])
        command_rows.append(
            {
                "family_category": family,
                "model_sha256": model["model_sha256"],
                "hmm_path": hmm_path.as_posix(),
                "target_fasta_path": target_fasta_path.as_posix(),
                "target_fasta_sha256": _sha256(target_fasta_path),
                "domtblout_path": domtblout_path.as_posix(),
                "main_output_path": main_output_path.as_posix(),
                "control_count": str(len(family_rows)),
                "hard_challenge_count": str(hard_count),
                "boundary_observation_count": str(sum(row["hard_negative"] == "no" for row in family_rows)),
                "command": _hmmsearch_command(hmm_path, target_fasta_path, domtblout_path, main_output_path),
                "command_status": COMMAND_STATUS,
                "notes": (
                    "Prepared from checksum-locked P05 models and controls; hmmsearch has not run. "
                    "Boundary observations are report-only and cannot set rejection thresholds."
                ),
            }
        )
    manifest_path = _write_tsv(
        calibration_dir / CALIBRATION_COMMAND_MANIFEST_FILENAME,
        command_rows,
        CALIBRATION_COMMAND_FIELDNAMES,
    )
    summary_rows = [
        {"kind": "total", "name": "calibration_jobs", "count": str(len(command_rows))},
        {"kind": "total", "name": "control_records", "count": str(sum(int(row["control_count"]) for row in command_rows))},
        {"kind": "total", "name": "hard_challenges", "count": str(sum(int(row["hard_challenge_count"]) for row in command_rows))},
        {"kind": "total", "name": "boundary_observations", "count": str(sum(int(row["boundary_observation_count"]) for row in command_rows))},
        {"kind": "command_status", "name": COMMAND_STATUS, "count": str(len(command_rows))},
    ]
    summary_path = _write_tsv(
        calibration_dir / CALIBRATION_COMMAND_SUMMARY_FILENAME,
        summary_rows,
        SUMMARY_FIELDNAMES,
    )
    return {"manifest": manifest_path, "summary": summary_path}


LEAVE_ONE_OUT_COMMAND_FIELDNAMES = (
    "family_category",
    "parent_model_sha256",
    "holdout_seed_id",
    "holdout_accession",
    "training_seed_count",
    "training_bundle_path",
    "training_bundle_sha256",
    "alignment_path",
    "leave_one_out_hmm_path",
    "positive_fasta_path",
    "positive_residue_sha256",
    "domtblout_path",
    "main_output_path",
    "command",
    "command_status",
    "notes",
)
LEAVE_ONE_OUT_RESULT_FIELDNAMES = (
    "family_category",
    "holdout_accession",
    "domtblout_path",
    "positive_hit_status",
    "best_full_score",
    "hmm_coverage",
    "domain_count",
)
CALIBRATION_DECISION_FIELDNAMES = (
    "family_category",
    "leave_one_out_variants",
    "positive_recovered",
    "positive_recovery_missing",
    "minimum_positive_full_score",
    "minimum_positive_hmm_coverage",
    "proposed_score_threshold",
    "proposed_hmm_coverage_threshold",
    "hard_challenge_count",
    "hard_challenge_hits",
    "hard_challenges_passing_proposed_rule",
    "boundary_observation_count",
    "boundary_observation_hits",
    "recommendation",
    "notes",
)
CONTROL_SMOKE_RESULT_FIELDNAMES = (
    "family_category",
    "control_id",
    "control_role",
    "hit_status",
    "best_full_score",
    "hmm_coverage",
    "domain_count",
)


def _load_seed_records_for_leave_one_out(path: Path, models: dict[str, dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    required = ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path")
    records_by_family: dict[str, list[dict[str, str]]] = {family: [] for family in models}
    seen_keys: set[tuple[str, str]] = set()
    for line_number, row in enumerate(_load_tsv(path, required), start=2):
        family = row["family_category"]
        if family not in models:
            continue
        accession = row["source_accession"]
        key = (family, accession)
        if not all(key) or not row["seed_id"]:
            raise ValueError(f"{path}:{line_number} has an empty seed identity")
        if key in seen_keys:
            raise ValueError(f"{path}:{line_number} duplicates profile seed {key!r}")
        seen_keys.add(key)
        if row["model_sha256"].lower() != models[family]["model_sha256"]:
            raise ValueError(f"{path}:{line_number} model SHA256 does not match the registry for {family}")
        sequence_path = Path(row["sequence_path"])
        if not sequence_path.is_file():
            raise FileNotFoundError(f"{path}:{line_number} profile seed FASTA is missing: {sequence_path}")
        normalized = dict(row)
        normalized["sequence_path"] = sequence_path.as_posix()
        records_by_family[family].append(normalized)
    for family, records in records_by_family.items():
        if len(records) < MINIMUM_LEAVE_ONE_OUT_SEEDS:
            raise ValueError(
                f"{family} has {len(records)} current profile seeds; at least {MINIMUM_LEAVE_ONE_OUT_SEEDS} are required "
                "for leave-one-out calibration with three retained training sequences"
            )
    return {
        family: sorted(records, key=lambda row: (row["seed_id"], row["source_accession"]))
        for family, records in records_by_family.items()
    }


def _write_seed_bundle(path: Path, records: list[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for record in records:
            handle.write(f">{record['seed_id']}|{record['source_accession']}\n")
            sequence = _read_single_fasta_sequence(Path(record["sequence_path"]))
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")
    return path


def _write_holdout_positive(path: Path, family: str, record: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    sequence = _read_single_fasta_sequence(Path(record["sequence_path"]))
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(f">positive|{family}|{record['source_accession']}\n")
        for start in range(0, len(sequence), 80):
            handle.write(f"{sequence[start:start + 80]}\n")
    return path


def build_leave_one_out_command_manifest(
    seed_registry: Path,
    model_registry: Path,
    calibration_dir: Path,
) -> dict[str, Path]:
    """Write deterministic leave-one-out MAFFT, hmmbuild, and positive-search commands."""

    models = _load_models_with_paths(model_registry)
    records_by_family = _load_seed_records_for_leave_one_out(seed_registry, models)
    root = calibration_dir / "leave_one_out"
    command_rows: list[dict[str, str]] = []
    for family, records in sorted(records_by_family.items()):
        for holdout in sorted(records, key=lambda row: row["source_accession"]):
            accession = holdout["source_accession"]
            variant_dir = root / family / accession
            training_records = [record for record in records if record["source_accession"] != accession]
            training_bundle_path = _write_seed_bundle(variant_dir / "training.faa", training_records)
            alignment_path = variant_dir / "training.aligned.faa"
            hmm_path = variant_dir / "leave_one_out.hmm"
            positive_fasta_path = _write_holdout_positive(variant_dir / "holdout_positive.faa", family, holdout)
            raw_dir = variant_dir / "raw_domtblout"
            log_dir = variant_dir / "hmmer_logs"
            raw_dir.mkdir(parents=True, exist_ok=True)
            log_dir.mkdir(parents=True, exist_ok=True)
            domtblout_path = raw_dir / "holdout_positive.domtblout"
            main_output_path = log_dir / "holdout_positive.txt"
            command_rows.append(
                {
                    "family_category": family,
                    "parent_model_sha256": models[family]["model_sha256"],
                    "holdout_seed_id": holdout["seed_id"],
                    "holdout_accession": accession,
                    "training_seed_count": str(len(training_records)),
                    "training_bundle_path": training_bundle_path.as_posix(),
                    "training_bundle_sha256": _sha256(training_bundle_path),
                    "alignment_path": alignment_path.as_posix(),
                    "leave_one_out_hmm_path": hmm_path.as_posix(),
                    "positive_fasta_path": positive_fasta_path.as_posix(),
                    "positive_residue_sha256": _sequence_residue_sha256(positive_fasta_path),
                    "domtblout_path": domtblout_path.as_posix(),
                    "main_output_path": main_output_path.as_posix(),
                    "command": (
                        "mafft --localpair --maxiterate 1000 --inputorder "
                        f"{_shell_quote(training_bundle_path.as_posix())} > {_shell_quote(alignment_path.as_posix())} && "
                        f"hmmbuild --amino {_shell_quote(hmm_path.as_posix())} {_shell_quote(alignment_path.as_posix())} && "
                        + _hmmsearch_command(hmm_path, positive_fasta_path, domtblout_path, main_output_path)
                    ),
                    "command_status": COMMAND_STATUS,
                    "notes": (
                        "Held-out seed is absent from the MAFFT/hmmbuild training set and is searched only after model build. "
                        "This is sequence-validation evidence, not a phenotype assay."
                    ),
                }
            )
    manifest_path = _write_tsv(
        calibration_dir / LEAVE_ONE_OUT_COMMAND_MANIFEST_FILENAME,
        command_rows,
        LEAVE_ONE_OUT_COMMAND_FIELDNAMES,
    )
    summary_rows = [
        {"kind": "total", "name": "leave_one_out_jobs", "count": str(len(command_rows))},
        {"kind": "total", "name": "families", "count": str(len(records_by_family))},
        {"kind": "command_status", "name": COMMAND_STATUS, "count": str(len(command_rows))},
    ]
    summary_path = _write_tsv(
        calibration_dir / LEAVE_ONE_OUT_COMMAND_SUMMARY_FILENAME,
        summary_rows,
        SUMMARY_FIELDNAMES,
    )
    return {"manifest": manifest_path, "summary": summary_path}


def _hmm_coverage(intervals: list[tuple[int, int]], hmm_length: int) -> float:
    """Return the union HMM-coordinate coverage for a reported target."""

    covered = 0
    previous_end = -1
    for start, end in sorted(intervals):
        if start < 1 or end < start or end > hmm_length:
            raise ValueError(f"Invalid HMM interval {start}-{end} for HMM length {hmm_length}")
        if start > previous_end:
            covered += end - start + 1
            previous_end = end
        elif end > previous_end:
            covered += end - previous_end
            previous_end = end
    return covered / hmm_length


def parse_leave_one_out_results(command_manifest: Path) -> list[dict[str, str]]:
    """Parse one held-out-positive HMMER result per checksum-locked variant.

    Each variant searches a single held-out FASTA record. A missing domtblout
    data line is therefore retained as an explicit failed positive recovery,
    rather than being silently omitted from calibration evidence.
    """

    manifest_rows = _load_tsv(
        command_manifest,
        ("family_category", "holdout_accession", "domtblout_path"),
    )
    results: list[dict[str, str]] = []
    seen_variants: set[tuple[str, str]] = set()
    for line_number, row in enumerate(manifest_rows, start=2):
        family = row["family_category"]
        accession = row["holdout_accession"]
        domtblout_path = Path(row["domtblout_path"])
        variant = (family, accession)
        if not all(variant):
            raise ValueError(f"{command_manifest}:{line_number} has an empty leave-one-out identity")
        if variant in seen_variants:
            raise ValueError(f"{command_manifest}:{line_number} duplicates leave-one-out variant {variant!r}")
        seen_variants.add(variant)
        if not domtblout_path.is_file():
            raise FileNotFoundError(f"{command_manifest}:{line_number} is missing domtblout: {domtblout_path}")

        expected_target = f"positive|{family}|{accession}"
        full_scores: list[float] = []
        intervals: list[tuple[int, int]] = []
        hmm_lengths: set[int] = set()
        for raw_line in domtblout_path.read_text(encoding="ascii").splitlines():
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.split()
            if len(fields) < 22:
                raise ValueError(f"{domtblout_path} has a malformed HMMER domtblout row")
            if fields[0] != expected_target:
                raise ValueError(
                    f"{domtblout_path} reports unexpected target {fields[0]!r}; expected {expected_target!r}"
                )
            try:
                hmm_length = int(fields[5])
                full_score = float(fields[7])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
            except ValueError as error:
                raise ValueError(f"{domtblout_path} has non-numeric HMMER fields") from error
            if hmm_length < 1:
                raise ValueError(f"{domtblout_path} has a nonpositive HMM length")
            hmm_lengths.add(hmm_length)
            full_scores.append(full_score)
            intervals.append((hmm_from, hmm_to))

        result = {
            "family_category": family,
            "holdout_accession": accession,
            "domtblout_path": domtblout_path.as_posix(),
            "positive_hit_status": "missing",
            "best_full_score": "",
            "hmm_coverage": "",
            "domain_count": "0",
        }
        if full_scores:
            if len(hmm_lengths) != 1:
                raise ValueError(f"{domtblout_path} reports inconsistent HMM lengths for one held-out positive")
            result.update(
                {
                    "positive_hit_status": "recovered",
                    "best_full_score": f"{max(full_scores):.1f}",
                    "hmm_coverage": f"{_hmm_coverage(intervals, hmm_lengths.pop()):.6f}",
                    "domain_count": str(len(full_scores)),
                }
            )
        results.append(result)
    return results


def derive_calibration_decisions(
    leave_one_out_results: list[dict[str, str]],
    control_results: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Derive conservative per-family recommendations from locked evidence.

    A proposed rule is the strictest score-and-coverage conjunction that still
    retains every leave-one-out positive. It is only an auditable candidate
    for human review and never changes a P06 approval field.
    """

    families = sorted({row["family_category"] for row in leave_one_out_results})
    decisions: list[dict[str, str]] = []
    for family in families:
        positives = [row for row in leave_one_out_results if row["family_category"] == family]
        recovered = [row for row in positives if row["positive_hit_status"] == "recovered"]
        missing_count = len(positives) - len(recovered)
        controls = [row for row in control_results if row["family_category"] == family]
        hard_controls = [row for row in controls if row["control_role"] == "cross_family_challenge"]
        boundary_controls = [row for row in controls if row["control_role"] == "boundary_observation"]
        if not hard_controls:
            raise ValueError(f"Calibration evidence lacks a hard cross-family challenge for {family}")
        hard_hits = [row for row in hard_controls if row["hit_status"] == "hit"]
        boundary_hits = [row for row in boundary_controls if row["hit_status"] == "hit"]
        decision = {
            "family_category": family,
            "leave_one_out_variants": str(len(positives)),
            "positive_recovered": str(len(recovered)),
            "positive_recovery_missing": str(missing_count),
            "minimum_positive_full_score": "",
            "minimum_positive_hmm_coverage": "",
            "proposed_score_threshold": "",
            "proposed_hmm_coverage_threshold": "",
            "hard_challenge_count": str(len(hard_controls)),
            "hard_challenge_hits": str(len(hard_hits)),
            "hard_challenges_passing_proposed_rule": "",
            "boundary_observation_count": str(len(boundary_controls)),
            "boundary_observation_hits": str(len(boundary_hits)),
            "recommendation": "blocked_positive_recovery_failed",
            "notes": "At least one held-out positive has no reportable HMMER domtblout match.",
        }
        if missing_count:
            decisions.append(decision)
            continue

        minimum_score = min(float(row["best_full_score"]) for row in recovered)
        minimum_coverage = min(float(row["hmm_coverage"]) for row in recovered)
        hard_passing = [
            row
            for row in hard_hits
            if float(row["best_full_score"]) >= minimum_score and float(row["hmm_coverage"]) >= minimum_coverage
        ]
        decision.update(
            {
                "minimum_positive_full_score": f"{minimum_score:.1f}",
                "minimum_positive_hmm_coverage": f"{minimum_coverage:.6f}",
                "proposed_score_threshold": f"{minimum_score:.1f}",
                "proposed_hmm_coverage_threshold": f"{minimum_coverage:.6f}",
                "hard_challenges_passing_proposed_rule": str(len(hard_passing)),
            }
        )
        if hard_passing:
            decision.update(
                {
                    "recommendation": "blocked_cross_family_overlap",
                    "notes": (
                        "One or more cross-family hard challenges pass the strictest score-and-coverage rule "
                        "that retains every leave-one-out positive."
                    ),
                }
            )
        else:
            decision.update(
                {
                    "recommendation": "eligible_for_human_review",
                    "notes": (
                        "All leave-one-out positives are retained and no hard cross-family challenge passes the "
                        "proposed rule; P06 remains blocked pending human review and registry update."
                    ),
                }
            )
        decisions.append(decision)
    return decisions


def parse_control_smoke_results(
    control_panel: Path,
    command_manifest: Path,
) -> list[dict[str, str]]:
    """Parse full-model controls while retaining no-hit panel records."""

    panel_rows = _load_tsv(control_panel, ("family_category", "control_id", "control_role"))
    command_rows = _load_tsv(command_manifest, ("family_category", "domtblout_path"))
    commands_by_family: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(command_rows, start=2):
        family = row["family_category"]
        if not family or family in commands_by_family:
            raise ValueError(f"{command_manifest}:{line_number} has an empty or duplicate model family")
        commands_by_family[family] = row

    panels_by_family: dict[str, list[dict[str, str]]] = {}
    for line_number, row in enumerate(panel_rows, start=2):
        family = row["family_category"]
        control_id = row["control_id"]
        if not family or not control_id:
            raise ValueError(f"{control_panel}:{line_number} has an empty control identity")
        if family not in commands_by_family:
            raise ValueError(f"{control_panel}:{line_number} has no matching HMMER command for {family}")
        panels_by_family.setdefault(family, []).append(row)

    results: list[dict[str, str]] = []
    for family, family_panel_rows in sorted(panels_by_family.items()):
        domtblout_path = Path(commands_by_family[family]["domtblout_path"])
        if not domtblout_path.is_file():
            raise FileNotFoundError(f"{command_manifest} is missing control domtblout for {family}: {domtblout_path}")
        valid_ids = {row["control_id"] for row in family_panel_rows}
        hits: dict[str, dict[str, object]] = {}
        for raw_line in domtblout_path.read_text(encoding="ascii").splitlines():
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.split()
            if len(fields) < 22:
                raise ValueError(f"{domtblout_path} has a malformed HMMER domtblout row")
            control_id = fields[0]
            if control_id not in valid_ids:
                raise ValueError(f"{domtblout_path} reports target absent from the control panel: {control_id!r}")
            try:
                hmm_length = int(fields[5])
                full_score = float(fields[7])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
            except ValueError as error:
                raise ValueError(f"{domtblout_path} has non-numeric HMMER fields") from error
            record = hits.setdefault(control_id, {"scores": [], "intervals": [], "hmm_lengths": set()})
            record["scores"].append(full_score)
            record["intervals"].append((hmm_from, hmm_to))
            record["hmm_lengths"].add(hmm_length)

        for panel_row in family_panel_rows:
            control_id = panel_row["control_id"]
            result = {
                "family_category": family,
                "control_id": control_id,
                "control_role": panel_row["control_role"],
                "hit_status": "no_hit",
                "best_full_score": "",
                "hmm_coverage": "",
                "domain_count": "0",
            }
            if control_id in hits:
                record = hits[control_id]
                hmm_lengths = record["hmm_lengths"]
                if len(hmm_lengths) != 1:
                    raise ValueError(f"{domtblout_path} reports inconsistent HMM lengths for {control_id}")
                result.update(
                    {
                        "hit_status": "hit",
                        "best_full_score": f"{max(record['scores']):.1f}",
                        "hmm_coverage": f"{_hmm_coverage(record['intervals'], next(iter(hmm_lengths))):.6f}",
                        "domain_count": str(len(record["scores"])),
                    }
                )
            results.append(result)
    return results


def write_calibration_result_tables(
    manifest_dir: Path,
    leave_one_out_results: list[dict[str, str]],
    control_results: list[dict[str, str]],
) -> dict[str, Path]:
    """Write compact, Git-trackable calibration evidence and recommendations."""

    decisions = derive_calibration_decisions(leave_one_out_results, control_results)
    return {
        "leave_one_out": _write_tsv(
            manifest_dir / LEAVE_ONE_OUT_RESULT_FILENAME,
            leave_one_out_results,
            LEAVE_ONE_OUT_RESULT_FIELDNAMES,
        ),
        "control_smoke": _write_tsv(
            manifest_dir / CONTROL_SMOKE_RESULT_FILENAME,
            control_results,
            CONTROL_SMOKE_RESULT_FIELDNAMES,
        ),
        "decisions": _write_tsv(
            manifest_dir / CALIBRATION_DECISION_FILENAME,
            decisions,
            CALIBRATION_DECISION_FIELDNAMES,
        ),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create the compact P05 HMM calibration control panel.")
    parser.add_argument(
        "--reference-manifest",
        type=Path,
        default=Path("01_reference_library/manifests/reference_library.seed_manifest.tsv"),
    )
    parser.add_argument(
        "--seed-registry",
        type=Path,
        default=Path("04_family_profiles/manifests/p05_hmm_seed_registry.tsv"),
    )
    parser.add_argument(
        "--model-registry",
        type=Path,
        default=Path("04_family_profiles/manifests/p05_hmm_model_registry.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("04_family_profiles/manifests") / CONTROL_PANEL_FILENAME,
    )
    parser.add_argument(
        "--build-commands",
        action="store_true",
        help="Materialize ignored calibration target FASTAs and a planned hmmsearch command manifest.",
    )
    parser.add_argument(
        "--build-leave-one-out",
        action="store_true",
        help="Materialize ignored leave-one-out training bundles and planned MAFFT/hmmbuild/hmmsearch commands.",
    )
    parser.add_argument(
        "--calibration-dir",
        type=Path,
        default=Path("04_family_profiles/calibration"),
        help="Ignored directory for calibration FASTAs, raw HMMER outputs, logs, and command manifests.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    rows = build_control_panel(args.reference_manifest, args.seed_registry, args.model_registry)
    output = write_control_panel(args.out, rows)
    print(f"Calibration control panel written: {output}")
    print(f"Controls: {len(rows)}")
    if args.build_commands:
        outputs = build_calibration_command_manifest(output, args.model_registry, args.calibration_dir)
        print(f"Calibration command manifest written: {outputs['manifest']}")
        print(f"Calibration command summary written: {outputs['summary']}")
    if args.build_leave_one_out:
        outputs = build_leave_one_out_command_manifest(args.seed_registry, args.model_registry, args.calibration_dir)
        print(f"Leave-one-out command manifest written: {outputs['manifest']}")
        print(f"Leave-one-out command summary written: {outputs['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
