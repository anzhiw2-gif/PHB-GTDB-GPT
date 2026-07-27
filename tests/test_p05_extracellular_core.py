from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest

from scripts import p05_extracellular_core as extracellular_core
from scripts import p05_hmm_calibration as calibration


class P05ExtracellularCoreTest(unittest.TestCase):
    def _write_tsv(self, path: Path, fields: tuple[str, ...] | list[str], rows: list[dict[str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _write_sequence(self, root: Path, accession: str) -> Path:
        path = root / "sequences" / f"{accession}.faa"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f">{accession}\nMSTNPKPQRIT{accession[-1]}\n", encoding="ascii")
        return path

    def _reference_row(self, family: str, accession: str, sequence_path: Path) -> dict[str, str]:
        return {
            "family_category": family,
            "source_accession": accession,
            "organism": f"Organism {accession}",
            "taxonomic_domain": "Bacteria",
            "evidence_level": "E1",
            "profile_seed_status": "approved",
            "sequence_path": sequence_path.as_posix(),
            "source_database": "TestDB",
            "source_release": "test-release",
            "source_version": "test-version",
            "retrieval_date": "2026-07-27",
            "source_url": f"https://example.test/{accession}",
            "doi": "10.0000/test",
            "pmid": "12345678",
            "pmcid": "",
            "literature_support_scope": "experimentally characterized test hydrolase",
            "notes": "test provenance",
        }

    def test_prepare_calibration_pools_extracellular_seeds_and_adds_close_hydrolase_challenge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            external = [
                ("extracellular_mcl_pha_dep", "E1"),
                ("extracellular_mcl_pha_dep", "E2"),
                ("extracellular_scl_pha_dep_type_I", "E3"),
                ("extracellular_scl_pha_dep_type_II", "E4"),
            ]
            internal = ("intracellular_mcl_pha_dep", "I1")
            reference_rows = [self._reference_row(family, accession, self._write_sequence(root, accession)) for family, accession in external]
            reference_rows.append(self._reference_row(*internal, self._write_sequence(root, internal[1])))
            reference_manifest = root / "reference.tsv"
            self._write_tsv(reference_manifest, calibration.REFERENCE_REQUIRED_FIELDS, reference_rows)

            seed_registry = root / "seed_registry.tsv"
            self._write_tsv(
                seed_registry,
                ["family_category", "source_accession", "seed_id", "sequence_path", "model_sha256"],
                [
                    {
                        "family_category": family,
                        "source_accession": accession,
                        "seed_id": f"seed-{accession.lower()}",
                        "sequence_path": next(row for row in reference_rows if row["source_accession"] == accession)["sequence_path"],
                        "model_sha256": "a" * 64,
                    }
                    for family, accession in [*external, internal]
                ],
            )

            close_control = self._reference_row("close_non_target_hydrolase", "C1", self._write_sequence(root, "C1"))
            controls_manifest = root / "controls.tsv"
            self._write_tsv(controls_manifest, calibration.REFERENCE_REQUIRED_FIELDS, [close_control])
            model_path = root / "core.hmm"
            model_path.write_text("HMMER3/f dummy\n", encoding="ascii")

            outputs = extracellular_core.prepare_core_calibration(
                seed_registry,
                reference_manifest,
                controls_manifest,
                model_path,
                root / "out",
            )

            with outputs["seed_registry"].open("r", encoding="utf-8", newline="") as handle:
                seed_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["source_accession"] for row in seed_rows], ["E1", "E2", "E3", "E4"])
            self.assertTrue(all(row["family_category"] == extracellular_core.CORE_FAMILY for row in seed_rows))

            with outputs["control_panel"].open("r", encoding="utf-8", newline="") as handle:
                controls = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {(row["source_accession"], row["control_role"]) for row in controls},
                {("I1", "cross_family_challenge"), ("C1", "close_non_target_hydrolase")},
            )
            self.assertTrue(all(row["hard_negative"] == "yes" for row in controls))

            with outputs["leave_one_out_manifest"].open("r", encoding="utf-8", newline="") as handle:
                leave_one_out = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(leave_one_out), 4)
            self.assertTrue(all(row["training_seed_count"] == "3" for row in leave_one_out))

    def test_command_line_writes_core_seed_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for family, accession in (
                ("extracellular_mcl_pha_dep", "E1"),
                ("extracellular_mcl_pha_dep", "E2"),
                ("extracellular_scl_pha_dep_type_I", "E3"),
                ("extracellular_scl_pha_dep_type_II", "E4"),
            ):
                sequence_path = self._write_sequence(root, accession)
                rows.append(
                    {
                        "family_category": family,
                        "source_accession": accession,
                        "sequence_path": sequence_path.as_posix(),
                    }
                )
            seed_registry = root / "seed_registry.tsv"
            self._write_tsv(seed_registry, ["family_category", "source_accession", "sequence_path"], rows)
            bundle = root / "core.faa"
            script = Path(__file__).parents[1] / "scripts" / "p05_extracellular_core.py"

            completed = subprocess.run(
                [sys.executable, str(script), "--seed-registry", str(seed_registry), "--bundle-path", str(bundle)],
                check=True,
                capture_output=True,
                text=True,
                cwd=root,
            )

            self.assertTrue(bundle.is_file())
            self.assertIn("Build command:", completed.stdout)
            self.assertTrue((root / "04_family_profiles" / "alignments").is_dir())
            self.assertTrue((root / "04_family_profiles" / "hmms").is_dir())

    def test_finalize_core_calibration_writes_named_compact_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            external = [
                ("extracellular_mcl_pha_dep", "E1"),
                ("extracellular_mcl_pha_dep", "E2"),
                ("extracellular_scl_pha_dep_type_I", "E3"),
                ("extracellular_scl_pha_dep_type_II", "E4"),
            ]
            internal = ("intracellular_mcl_pha_dep", "I1")
            reference_rows = [self._reference_row(family, accession, self._write_sequence(root, accession)) for family, accession in external]
            reference_rows.append(self._reference_row(*internal, self._write_sequence(root, internal[1])))
            reference_manifest = root / "reference.tsv"
            self._write_tsv(reference_manifest, calibration.REFERENCE_REQUIRED_FIELDS, reference_rows)
            seed_registry = root / "seed_registry.tsv"
            self._write_tsv(
                seed_registry,
                ["family_category", "source_accession", "seed_id", "sequence_path", "model_sha256"],
                [
                    {
                        "family_category": family,
                        "source_accession": accession,
                        "seed_id": f"seed-{accession.lower()}",
                        "sequence_path": next(row for row in reference_rows if row["source_accession"] == accession)["sequence_path"],
                        "model_sha256": "a" * 64,
                    }
                    for family, accession in [*external, internal]
                ],
            )
            controls_manifest = root / "controls.tsv"
            self._write_tsv(
                controls_manifest,
                calibration.REFERENCE_REQUIRED_FIELDS,
                [self._reference_row("close_non_target_hydrolase", "C1", self._write_sequence(root, "C1"))],
            )
            model_path = root / "core.hmm"
            model_path.write_text("HMMER3/f dummy\n", encoding="ascii")
            outputs = extracellular_core.prepare_core_calibration(
                seed_registry, reference_manifest, controls_manifest, model_path, root / "out"
            )
            with outputs["leave_one_out_manifest"].open("r", encoding="utf-8", newline="") as handle:
                leave_rows = list(csv.DictReader(handle, delimiter="\t"))
            for row in leave_rows:
                Path(row["domtblout_path"]).write_text(
                    f"positive|{extracellular_core.CORE_FAMILY}|{row['holdout_accession']} - 100 core - 100 1e-20 50.0 0.0 1 1 1e-20 1e-20 50.0 0.0 1 100 1 100 1 100 0.99 fixture\n",
                    encoding="ascii",
                )
            with outputs["control_manifest"].open("r", encoding="utf-8", newline="") as handle:
                control_manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            for row in control_manifest_rows:
                Path(row["domtblout_path"]).write_text("# no controls pass\n", encoding="ascii")

            compact = extracellular_core.finalize_core_calibration(root / "out", root / "manifests")

            with compact["decisions"].open("r", encoding="utf-8", newline="") as handle:
                decisions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(decisions[0]["positive_recovered"], "4")
            self.assertEqual(decisions[0]["hard_challenge_count"], "2")
            self.assertEqual(decisions[0]["recommendation"], "eligible_for_human_review")


if __name__ == "__main__":
    unittest.main()
