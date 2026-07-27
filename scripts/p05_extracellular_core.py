"""Prepare a pooled extracellular PHA-depolymerase HMM calibration package.

The three experimentally supported extracellular PHA-depolymerase groups share
a broad alpha/beta-hydrolase core and failed mutually exclusive HMM calibration.
This helper pools their verified seeds for a conservative top-level scan model.
Its subtype labels remain provenance annotations for downstream architecture and
phylogenetic review; they are not inferred from the pooled HMM hit alone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

try:
    from scripts import p05_hmm_calibration as calibration
except ModuleNotFoundError:  # Supports ``python scripts/p05_extracellular_core.py`` on T141.
    import p05_hmm_calibration as calibration


CORE_FAMILY = "extracellular_pha_depolymerase_core"
EXTRACELLULAR_SUBFAMILIES = frozenset(
    {
        "extracellular_mcl_pha_dep",
        "extracellular_scl_pha_dep_type_I",
        "extracellular_scl_pha_dep_type_II",
    }
)
CORE_SEED_FIELDS = ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path")
CORE_MODEL_FIELDS = ("family_category", "model_sha256", "model_path")


def _read_tsv(path: Path, required_fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing = [field for field in required_fields if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def _write_tsv(path: Path, rows: list[dict[str, str]], fields: tuple[str, ...]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _core_seed_records(seed_registry: Path) -> list[dict[str, str]]:
    rows = _read_tsv(seed_registry, ("family_category", "source_accession", "sequence_path"))
    selected = [row for row in rows if row["family_category"] in EXTRACELLULAR_SUBFAMILIES]
    if len(selected) < calibration.MINIMUM_LEAVE_ONE_OUT_SEEDS:
        raise ValueError("The extracellular core needs at least four current experimentally supported seeds")
    seen: set[str] = set()
    for row in selected:
        accession = row["source_accession"]
        if not accession or accession in seen:
            raise ValueError(f"The extracellular core seed registry has a duplicate or empty accession: {accession!r}")
        if not Path(row["sequence_path"]).is_file():
            raise FileNotFoundError(f"Extracellular core seed FASTA is missing: {row['sequence_path']}")
        seen.add(accession)
    return sorted(selected, key=lambda row: row["source_accession"])


def write_core_seed_bundle(seed_registry: Path, bundle_path: Path) -> Path:
    """Materialize the deterministic pooled seed bundle before HMM construction."""

    records = _core_seed_records(seed_registry)
    bundle_path.parent.mkdir(parents=True, exist_ok=True)
    with bundle_path.open("w", encoding="ascii", newline="\n") as handle:
        for record in records:
            accession = record["source_accession"]
            handle.write(f">core_seed|{record['family_category']}|{accession}\n")
            sequence = calibration._read_single_fasta_sequence(Path(record["sequence_path"]))
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")
    return bundle_path


def _reference_index(reference_manifest: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = _read_tsv(reference_manifest, calibration.REFERENCE_REQUIRED_FIELDS)
    index = {(row["family_category"], row["source_accession"]): row for row in rows}
    if len(index) != len(rows):
        raise ValueError(f"{reference_manifest} contains duplicate family/accession records")
    return index


def _control_row(
    source: dict[str, str],
    *,
    model_sha256: str,
    role: str,
    control_family: str,
) -> dict[str, str]:
    accession = source["source_accession"]
    sequence_path = Path(source["sequence_path"])
    return {
        "family_category": CORE_FAMILY,
        "model_sha256": model_sha256,
        "control_id": f"{CORE_FAMILY}|{role}|{control_family}|{accession}",
        "control_role": role,
        "hard_negative": "yes",
        "expected_outcome": "must_fail_threshold",
        "control_family_category": control_family,
        "source_accession": accession,
        "organism": source["organism"],
        "taxonomic_domain": source["taxonomic_domain"],
        "evidence_level": source["evidence_level"],
        "profile_seed_status": source["profile_seed_status"],
        "sequence_path": sequence_path.as_posix(),
        "sequence_sha256": calibration._sequence_residue_sha256(sequence_path),
        "source_database": source["source_database"],
        "source_release": source["source_release"],
        "source_version": source["source_version"],
        "retrieval_date": source["retrieval_date"],
        "source_url": source["source_url"],
        "doi": source["doi"],
        "pmid": source["pmid"],
        "pmcid": source["pmcid"],
        "literature_support_scope": source["literature_support_scope"],
        "architecture_rationale": (
            "Hard specificity challenge for the pooled extracellular PHA-depolymerase core; "
            "a passing final rule blocks P06 approval."
        ),
        "notes": source["notes"],
    }


def prepare_core_calibration(
    seed_registry: Path,
    reference_manifest: Path,
    controls_manifest: Path,
    model_path: Path,
    output_dir: Path,
) -> dict[str, Path]:
    """Create checksum-locked core inputs compatible with the generic P05 parser."""

    if not model_path.is_file():
        raise FileNotFoundError(f"Pooled extracellular core HMM is missing: {model_path}")
    model_sha256 = _sha256(model_path)
    core_seeds = _core_seed_records(seed_registry)
    reference_by_key = _reference_index(reference_manifest)
    all_seed_rows = _read_tsv(seed_registry, ("family_category", "source_accession", "sequence_path"))
    close_controls = _read_tsv(controls_manifest, calibration.REFERENCE_REQUIRED_FIELDS)

    output_dir.mkdir(parents=True, exist_ok=True)
    core_seed_rows = [
        {
            "family_category": CORE_FAMILY,
            "model_sha256": model_sha256,
            "seed_id": f"core-{record['source_accession']}",
            "source_accession": record["source_accession"],
            "sequence_path": record["sequence_path"],
        }
        for record in core_seeds
    ]
    core_seed_registry = _write_tsv(output_dir / "core_seed_registry.tsv", core_seed_rows, CORE_SEED_FIELDS)
    core_model_registry = _write_tsv(
        output_dir / "core_model_registry.tsv",
        [{"family_category": CORE_FAMILY, "model_sha256": model_sha256, "model_path": model_path.as_posix()}],
        CORE_MODEL_FIELDS,
    )

    control_rows: list[dict[str, str]] = []
    seen_control_accessions: set[str] = set()
    for record in sorted(all_seed_rows, key=lambda row: (row["family_category"], row["source_accession"])):
        family = record["family_category"]
        if family in EXTRACELLULAR_SUBFAMILIES:
            continue
        key = (family, record["source_accession"])
        source = reference_by_key.get(key)
        if source is None:
            raise ValueError(f"No reference provenance for core cross-family challenge {key!r}")
        control_rows.append(
            _control_row(source, model_sha256=model_sha256, role="cross_family_challenge", control_family=family)
        )
        seen_control_accessions.add(source["source_accession"])
    for source in sorted(close_controls, key=lambda row: row["source_accession"]):
        accession = source["source_accession"]
        if not accession or accession in seen_control_accessions:
            raise ValueError(f"Close non-target control accession is empty or duplicated: {accession!r}")
        if not Path(source["sequence_path"]).is_file():
            raise FileNotFoundError(f"Close non-target control FASTA is missing: {source['sequence_path']}")
        control_rows.append(
            _control_row(
                source,
                model_sha256=model_sha256,
                role="close_non_target_hydrolase",
                control_family=source["family_category"],
            )
        )
        seen_control_accessions.add(accession)

    if not close_controls:
        raise ValueError("At least one experimentally characterized close non-target hydrolase is required")
    control_panel = calibration.write_control_panel(output_dir / "core_control_panel.tsv", control_rows)
    calibration_dir = output_dir / "calibration"
    control_outputs = calibration.build_calibration_command_manifest(control_panel, core_model_registry, calibration_dir)
    leave_one_out_outputs = calibration.build_leave_one_out_command_manifest(core_seed_registry, core_model_registry, calibration_dir)
    return {
        "seed_registry": core_seed_registry,
        "model_registry": core_model_registry,
        "control_panel": control_panel,
        "control_manifest": control_outputs["manifest"],
        "leave_one_out_manifest": leave_one_out_outputs["manifest"],
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare pooled extracellular PHA-depolymerase HMM inputs.")
    parser.add_argument("--seed-registry", type=Path, default=Path("04_family_profiles/manifests/p05_hmm_seed_registry.tsv"))
    parser.add_argument("--reference-manifest", type=Path, default=Path("01_reference_library/reference_library.normalized.tsv"))
    parser.add_argument("--controls-manifest", type=Path, default=Path("04_family_profiles/manifests/p05_extracellular_core_close_controls.tsv"))
    parser.add_argument("--bundle-path", type=Path, default=Path("04_family_profiles/seed_bundles/extracellular_pha_depolymerase_core.faa"))
    parser.add_argument("--alignment-path", type=Path, default=Path("04_family_profiles/alignments/extracellular_pha_depolymerase_core.aligned.faa"))
    parser.add_argument("--model-path", type=Path, default=Path("04_family_profiles/hmms/extracellular_pha_depolymerase_core.hmm"))
    parser.add_argument("--output-dir", type=Path, default=Path("04_family_profiles/calibration/extracellular_pha_depolymerase_core"))
    parser.add_argument("--prepare-calibration", action="store_true", help="Write checksum-locked control and leave-one-out manifests after hmmbuild.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.prepare_calibration:
        outputs = prepare_core_calibration(
            args.seed_registry,
            args.reference_manifest,
            args.controls_manifest,
            args.model_path,
            args.output_dir,
        )
        for name, path in sorted(outputs.items()):
            print(f"{name}: {path}")
        return 0

    bundle_path = write_core_seed_bundle(args.seed_registry, args.bundle_path)
    print(f"Seed bundle written: {bundle_path}")
    print(
        "Build command: "
        f"mafft --localpair --maxiterate 1000 --inputorder {calibration._shell_quote(bundle_path.as_posix())} "
        f"> {calibration._shell_quote(args.alignment_path.as_posix())} && "
        f"hmmbuild --amino {calibration._shell_quote(args.model_path.as_posix())} "
        f"{calibration._shell_quote(args.alignment_path.as_posix())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
