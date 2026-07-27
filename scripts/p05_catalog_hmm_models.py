"""Create compact, Git-trackable provenance for locally built P05 HMMs.

The HMMs, seed bundles, and alignments are deliberately machine-local. This
script records their deterministic SHA256 values, HMMER header metadata, and
the accession-level reference provenance needed to decide whether a profile is
eligible for a later P06 scan. It does not build profiles or make phenotype
claims from sequence similarity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


MODEL_REGISTRY_FILENAME = "p05_hmm_model_registry.tsv"
SEED_REGISTRY_FILENAME = "p05_hmm_seed_registry.tsv"
PROPOSED_UPDATES_FILENAME = "p05_hmm_proposed_seed_updates.tsv"

MODEL_FIELDNAMES = (
    "family_category",
    "model_status",
    "approved_for_p06",
    "scan_permission",
    "model_path",
    "model_sha256",
    "hmmer_format",
    "hmmer_version",
    "model_name",
    "model_length",
    "seed_sequence_count",
    "effective_sequence_number",
    "model_build_date",
    "model_specific_thresholds",
    "seed_bundle_path",
    "seed_bundle_sha256",
    "alignment_path",
    "alignment_sha256",
    "recorded_hmmbuild_command",
    "calibration_status",
    "seed_decision_document",
    "notes",
)
SEED_FIELDNAMES = (
    "family_category",
    "model_sha256",
    "current_bundle_role",
    "next_rebuild_decision",
    "seed_id",
    "source_accession",
    "organism",
    "taxonomic_domain",
    "evidence_level",
    "sequence_length_aa",
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
    "bundle_header",
    "bundle_sha256",
    "bundle_membership_verified",
    "notes",
)
PROPOSED_FIELDNAMES = (
    "family_category",
    "accession",
    "organism",
    "evidence_level",
    "proposed_role",
    "decision_status",
    "experimental_or_architecture_support",
    "doi",
    "pmid",
    "source_database",
    "retrieval_status",
    "sequence_sha256",
    "notes",
)

DECISION_DOCUMENT = "docs/P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md"

MODEL_STATUS = {
    "archaeal_patatin_like_pha_dep": "blocked_for_rebuild",
    "intracellular_mcl_pha_dep": "provisional_archived_needs_row_audit",
    "intracellular_phaZ_no_lipase_box": "proposed_seed_update_pending_user_confirmation",
}

SEED_DECISIONS = {
    "archaeal_patatin_like_pha_dep": {
        "AFK21580.1": "retain_experimental_patin_like_anchor",
        "CCQ36014.1": "retain_e3_patin_like_coverage",
        "AHB64615.1": "demote_to_non_patin_architecture_boundary",
        "AHZ23723.1": "demote_to_non_patin_architecture_boundary",
        "AJF25805.1": "demote_to_non_patin_architecture_boundary",
        "EFW93255.1": "demote_to_non_patin_architecture_boundary",
        "KOX95185.1": "demote_to_non_patin_architecture_boundary",
        "KZX50211.1": "demote_to_non_patin_architecture_boundary",
    },
    "intracellular_mcl_pha_dep": {
        "Q5Y152": "retain_primary_experimental_anchor",
    },
    "intracellular_phaZ_no_lipase_box": {
        "O87189": "retain_experimental_profile_seed",
        "Q0K7T2": "retain_experimental_profile_seed",
        "Q71KW6": "retain_experimental_profile_seed",
        "Q92TD3": "retain_experimental_profile_seed",
        "Q0K4D5": "demote_to_boundary_control",
    },
}

PROPOSED_UPDATES = (
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "AFK21580.1",
        "organism": "Haloferax mediterranei ATCC 33500",
        "evidence_level": "E1",
        "proposed_role": "retain_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "PhaZh1/HFX_6464 experimental anchor; patatin PF01734 and N-terminal GTSGG motif",
        "doi": "10.1128/AEM.04269-14",
        "pmid": "25710370",
        "source_database": "NCBI Protein",
        "retrieval_status": "already_retrieved_in_current_bundle",
        "sequence_sha256": "",
        "notes": "Retain only in a coherent PhaZh1-like patatin profile.",
    },
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "CCQ36014.1",
        "organism": "Natronomonas moolapensis 8.8.11",
        "evidence_level": "E3",
        "proposed_role": "retain_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Patatin/RssA-related architecture; 323 aa; N-terminal GTSGG motif",
        "doi": "",
        "pmid": "",
        "source_database": "NCBI Protein",
        "retrieval_status": "already_retrieved_in_current_bundle",
        "sequence_sha256": "",
        "notes": "Coverage evidence only, not phenotype evidence.",
    },
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "CCQ32286.1",
        "organism": "Halorhabdus tiamatea",
        "evidence_level": "E3",
        "proposed_role": "add_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Patatin/RssA-related candidate; 322 aa; N-terminal GSSGG motif",
        "doi": "",
        "pmid": "",
        "source_database": "NCBI Protein",
        "retrieval_status": "not_retrieved_pending_confirmation",
        "sequence_sha256": "",
        "notes": "Coverage evidence only, not phenotype evidence.",
    },
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "AGN01047.1",
        "organism": "Salinarchaeum sp. Harcht-Bsk1",
        "evidence_level": "E3",
        "proposed_role": "add_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Patatin/RssA-related candidate; 321 aa; N-terminal GTSGG motif",
        "doi": "",
        "pmid": "",
        "source_database": "NCBI Protein",
        "retrieval_status": "not_retrieved_pending_confirmation",
        "sequence_sha256": "",
        "notes": "Coverage evidence only, not phenotype evidence.",
    },
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "KYH27761.1",
        "organism": "Halalkalicoccus paucihalophilus",
        "evidence_level": "E3",
        "proposed_role": "add_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Patatin/RssA-related candidate; 329 aa; N-terminal GTSGG motif",
        "doi": "",
        "pmid": "",
        "source_database": "NCBI Protein",
        "retrieval_status": "not_retrieved_pending_confirmation",
        "sequence_sha256": "",
        "notes": "Coverage evidence only, not phenotype evidence.",
    },
    {
        "family_category": "archaeal_patatin_like_pha_dep",
        "accession": "ELY43313.1",
        "organism": "Natronorubrum tibetense",
        "evidence_level": "E3",
        "proposed_role": "add_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Patatin/RssA-related candidate; 269 aa; N-terminal GTSGG motif",
        "doi": "",
        "pmid": "",
        "source_database": "NCBI Protein",
        "retrieval_status": "not_retrieved_pending_confirmation",
        "sequence_sha256": "",
        "notes": "Coverage evidence only, not phenotype evidence.",
    },
    {
        "family_category": "intracellular_phaZ_no_lipase_box",
        "accession": "Q9WX79",
        "organism": "Paracoccus denitrificans",
        "evidence_level": "E2",
        "proposed_role": "add_profile_seed",
        "decision_status": "user_confirmation_required",
        "experimental_or_architecture_support": "Functional heterologous-expression PHB-granule degradation assay; PF06850/TIGR01849 family space",
        "doi": "10.1111/j.1574-6968.2001.tb10558.x",
        "pmid": "11267773",
        "source_database": "UniProtKB",
        "retrieval_status": "not_retrieved_pending_confirmation",
        "sequence_sha256": "",
        "notes": "Replace boundary paralog Q0K4D5 for cross-genus profile coverage.",
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_hmm_header(path: Path) -> dict[str, str]:
    """Parse the compact metadata fields emitted before an HMMER model body."""

    metadata: dict[str, str] = {"thresholds": "none"}
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped == "HMM":
                break
            if stripped.startswith("HMMER"):
                metadata["hmmer_format"] = stripped.split()[0]
                version_match = re.search(r"\[([^|\]]+)", stripped)
                metadata["hmmer_version"] = version_match.group(1).strip() if version_match else ""
                continue
            if not stripped or " " not in stripped:
                continue
            key, value = stripped.split(None, 1)
            if key in {"NAME", "LENG", "NSEQ", "EFFN", "DATE"}:
                metadata[key] = value.strip()
            if key in {"GA", "TC", "NC"}:
                metadata["thresholds"] = "present"
    required = ("hmmer_format", "hmmer_version", "NAME", "LENG", "NSEQ", "EFFN", "DATE")
    missing = [key for key in required if not metadata.get(key)]
    if missing:
        raise ValueError(f"{path} is missing HMMER header fields: {', '.join(missing)}")
    return metadata


def parse_bundle_headers(path: Path) -> list[tuple[str, str]]:
    headers: list[tuple[str, str]] = []
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            fields = header.split("|", maxsplit=1)
            if len(fields) != 2 or not all(fields):
                raise ValueError(f"{path} has an invalid seed bundle header: {header!r}")
            headers.append((fields[0], fields[1]))
    if not headers:
        raise ValueError(f"{path} does not contain FASTA headers")
    return headers


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def model_status(family: str) -> str:
    return MODEL_STATUS.get(family, "provisional_archived_pending_calibration")


def seed_decision(family: str, accession: str) -> str:
    explicit = SEED_DECISIONS.get(family, {}).get(accession)
    if explicit:
        return explicit
    if family == "intracellular_mcl_pha_dep":
        return "retain_only_after_accession_level_experimental_evidence_audit"
    return "current_provisional_profile_seed"


def catalog_hmm_models(
    reference_manifest: Path,
    command_manifest: Path,
    model_dir: Path,
    bundle_dir: Path,
    alignment_dir: Path,
    outdir: Path,
) -> dict[str, Path]:
    """Validate local generated artifacts and write compact tracked registries."""

    references = load_tsv(reference_manifest)
    required_reference_fields = {"family_category", "seed_id", "source_accession", "sequence_path"}
    missing_reference_fields = required_reference_fields - set(references[0]) if references else required_reference_fields
    if missing_reference_fields:
        raise ValueError(f"{reference_manifest} is missing required columns: {', '.join(sorted(missing_reference_fields))}")
    reference_by_key = {(row["family_category"], row["source_accession"]): row for row in references}
    commands = {row["family_category"]: row for row in load_tsv(command_manifest)}

    model_rows: list[dict[str, str]] = []
    seed_rows: list[dict[str, str]] = []
    for model_path in sorted(model_dir.glob("*.hmm")):
        family = model_path.stem
        bundle_path = bundle_dir / f"{family}.faa"
        alignment_path = alignment_dir / f"{family}.aligned.faa"
        if not bundle_path.is_file() or not alignment_path.is_file():
            raise ValueError(f"{family} is missing its seed bundle or alignment")
        header = parse_hmm_header(model_path)
        bundle_headers = parse_bundle_headers(bundle_path)
        if int(header["NSEQ"]) != len(bundle_headers):
            raise ValueError(f"{family} HMM NSEQ does not match its bundle header count")
        model_hash = sha256(model_path)
        bundle_hash = sha256(bundle_path)
        command = commands.get(family, {}).get("hmmbuild_command", "not_recorded")
        model_rows.append(
            {
                "family_category": family,
                "model_status": model_status(family),
                "approved_for_p06": "no",
                "scan_permission": "blocked",
                "model_path": model_path.as_posix(),
                "model_sha256": model_hash,
                "hmmer_format": header["hmmer_format"],
                "hmmer_version": header["hmmer_version"],
                "model_name": header["NAME"],
                "model_length": header["LENG"],
                "seed_sequence_count": header["NSEQ"],
                "effective_sequence_number": header["EFFN"],
                "model_build_date": header["DATE"],
                "model_specific_thresholds": header["thresholds"],
                "seed_bundle_path": bundle_path.as_posix(),
                "seed_bundle_sha256": bundle_hash,
                "alignment_path": alignment_path.as_posix(),
                "alignment_sha256": sha256(alignment_path),
                "recorded_hmmbuild_command": command,
                "calibration_status": "not_complete",
                "seed_decision_document": DECISION_DOCUMENT,
                "notes": "The raw model is machine-local and ignored by Git; this row is metadata only.",
            }
        )
        for seed_id, accession in bundle_headers:
            reference = reference_by_key.get((family, accession))
            if reference is None:
                raise ValueError(f"{family} bundle accession {accession!r} is absent from the reference manifest")
            sequence_path = Path(reference["sequence_path"])
            if not sequence_path.is_file():
                raise ValueError(f"{family} reference FASTA is missing: {sequence_path}")
            seed_rows.append(
                {
                    "family_category": family,
                    "model_sha256": model_hash,
                    "current_bundle_role": "provisional_model_seed",
                    "next_rebuild_decision": seed_decision(family, accession),
                    "seed_id": seed_id,
                    "source_accession": accession,
                    "organism": reference.get("organism", ""),
                    "taxonomic_domain": reference.get("taxonomic_domain", ""),
                    "evidence_level": reference.get("evidence_level", ""),
                    "sequence_length_aa": reference.get("sequence_length_aa", ""),
                    "sequence_path": reference["sequence_path"],
                    "sequence_sha256": sha256(sequence_path),
                    "source_database": reference.get("source_database", ""),
                    "source_release": reference.get("source_release", ""),
                    "source_version": reference.get("source_version", ""),
                    "retrieval_date": reference.get("retrieval_date", ""),
                    "source_url": reference.get("source_url", ""),
                    "doi": reference.get("doi", ""),
                    "pmid": reference.get("pmid", ""),
                    "pmcid": reference.get("pmcid", ""),
                    "literature_support_scope": reference.get("literature_support_scope", ""),
                    "bundle_header": f">{seed_id}|{accession}",
                    "bundle_sha256": bundle_hash,
                    "bundle_membership_verified": "yes",
                    "notes": reference.get("notes", ""),
                }
            )
    if not model_rows:
        raise ValueError(f"No HMM files were found in {model_dir}")
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "model_registry": outdir / MODEL_REGISTRY_FILENAME,
        "seed_registry": outdir / SEED_REGISTRY_FILENAME,
        "proposed_updates": outdir / PROPOSED_UPDATES_FILENAME,
    }
    write_tsv(outputs["model_registry"], model_rows, MODEL_FIELDNAMES)
    write_tsv(outputs["seed_registry"], seed_rows, SEED_FIELDNAMES)
    write_tsv(outputs["proposed_updates"], list(PROPOSED_UPDATES), PROPOSED_FIELDNAMES)
    return outputs


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write compact provenance registries for local P05 HMM artifacts.")
    parser.add_argument("--reference-manifest", type=Path, default=Path("01_reference_library/manifests/reference_library.seed_manifest.tsv"))
    parser.add_argument("--command-manifest", type=Path, default=Path("04_family_profiles/manifests/p05_family_profile_command_manifest.tsv"))
    parser.add_argument("--model-dir", type=Path, default=Path("04_family_profiles/hmms"))
    parser.add_argument("--bundle-dir", type=Path, default=Path("04_family_profiles/seed_bundles"))
    parser.add_argument("--alignment-dir", type=Path, default=Path("04_family_profiles/alignments"))
    parser.add_argument("--outdir", type=Path, default=Path("04_family_profiles/manifests"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    outputs = catalog_hmm_models(
        args.reference_manifest,
        args.command_manifest,
        args.model_dir,
        args.bundle_dir,
        args.alignment_dir,
        args.outdir,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
