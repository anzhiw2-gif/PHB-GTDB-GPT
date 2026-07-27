from __future__ import annotations

import csv
import hashlib
import tempfile
from pathlib import Path
import unittest

from scripts import p05_hmm_calibration as calibration


class P05HmmCalibrationTest(unittest.TestCase):
    def _write_tsv(self, path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _write_sequence(self, root: Path, accession: str) -> Path:
        path = root / f"{accession}.faa"
        path.write_text(f">seed|{accession}\nMSTNPKPQRIT\n", encoding="ascii")
        return path

    def _reference_row(
        self,
        *,
        family: str,
        accession: str,
        sequence_path: Path,
        profile_seed_status: str,
    ) -> dict[str, str]:
        return {
            "family_category": family,
            "source_accession": accession,
            "organism": f"Test organism {accession}",
            "taxonomic_domain": "Bacteria",
            "evidence_level": "E1",
            "profile_seed_status": profile_seed_status,
            "sequence_path": sequence_path.as_posix(),
            "source_database": "TestDB",
            "source_release": "test-release",
            "source_version": "test-version",
            "retrieval_date": "2026-07-27",
            "source_url": f"https://example.test/{accession}",
            "doi": "10.0000/test",
            "pmid": "12345678",
            "pmcid": "",
            "literature_support_scope": "test support scope",
            "notes": "test notes",
        }

    def _build_fixture(self, root: Path) -> tuple[Path, Path, Path]:
        target_seed = self._write_sequence(root, "A1")
        other_seed = self._write_sequence(root, "B1")
        boundary = self._write_sequence(root, "A2")
        reference_manifest = root / "reference.tsv"
        fields = list(calibration.REFERENCE_REQUIRED_FIELDS)
        self._write_tsv(
            reference_manifest,
            fields,
            [
                self._reference_row(
                    family="family_a",
                    accession="A1",
                    sequence_path=target_seed,
                    profile_seed_status="approved",
                ),
                self._reference_row(
                    family="family_a",
                    accession="A2",
                    sequence_path=boundary,
                    profile_seed_status="boundary_candidate",
                ),
                self._reference_row(
                    family="family_b",
                    accession="B1",
                    sequence_path=other_seed,
                    profile_seed_status="approved",
                ),
            ],
        )
        seed_registry = root / "seed_registry.tsv"
        self._write_tsv(
            seed_registry,
            ["family_category", "source_accession", "model_sha256"],
            [
                {"family_category": "family_a", "source_accession": "A1", "model_sha256": "a" * 64},
                {"family_category": "family_b", "source_accession": "B1", "model_sha256": "b" * 64},
            ],
        )
        model_registry = root / "model_registry.tsv"
        self._write_tsv(
            model_registry,
            ["family_category", "model_sha256"],
            [
                {"family_category": "family_a", "model_sha256": "a" * 64},
                {"family_category": "family_b", "model_sha256": "b" * 64},
            ],
        )
        return reference_manifest, seed_registry, model_registry

    def test_control_panel_uses_cross_family_seeds_as_hard_challenges_and_target_boundaries_as_observations(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            reference_manifest, seed_registry, model_registry = self._build_fixture(Path(tmp))

            rows = calibration.build_control_panel(reference_manifest, seed_registry, model_registry)

            by_key = {(row["family_category"], row["source_accession"]): row for row in rows}
            self.assertEqual(by_key[("family_a", "B1")]["control_role"], "cross_family_challenge")
            self.assertEqual(by_key[("family_a", "B1")]["hard_negative"], "yes")
            self.assertEqual(by_key[("family_a", "B1")]["expected_outcome"], "must_fail_threshold")
            self.assertEqual(by_key[("family_a", "A2")]["control_role"], "boundary_observation")
            self.assertEqual(by_key[("family_a", "A2")]["hard_negative"], "no")
            self.assertEqual(by_key[("family_a", "A2")]["expected_outcome"], "report_only")
            self.assertNotIn(("family_a", "A1"), by_key)
            self.assertEqual(len(by_key[("family_a", "B1")]["sequence_sha256"]), 64)

    def test_control_panel_rejects_a_model_without_cross_family_challenge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target_seed = self._write_sequence(root, "A1")
            reference_manifest = root / "reference.tsv"
            self._write_tsv(
                reference_manifest,
                list(calibration.REFERENCE_REQUIRED_FIELDS),
                [
                    self._reference_row(
                        family="family_a",
                        accession="A1",
                        sequence_path=target_seed,
                        profile_seed_status="approved",
                    )
                ],
            )
            seed_registry = root / "seed_registry.tsv"
            self._write_tsv(
                seed_registry,
                ["family_category", "source_accession", "model_sha256"],
                [{"family_category": "family_a", "source_accession": "A1", "model_sha256": "a" * 64}],
            )
            model_registry = root / "model_registry.tsv"
            self._write_tsv(
                model_registry,
                ["family_category", "model_sha256"],
                [{"family_category": "family_a", "model_sha256": "a" * 64}],
            )

            with self.assertRaisesRegex(ValueError, "cross-family challenge"):
                calibration.build_control_panel(reference_manifest, seed_registry, model_registry)

    def test_calibration_commands_materialize_checksum_locked_control_fastas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_manifest, seed_registry, model_registry = self._build_fixture(root)
            hmm_dir = root / "hmms"
            hmm_dir.mkdir()
            model_rows: list[dict[str, str]] = []
            for family in ("family_a", "family_b"):
                model_path = hmm_dir / f"{family}.hmm"
                model_path.write_text(f"HMMER3/f {family}\n", encoding="ascii")
                model_rows.append(
                    {
                        "family_category": family,
                        "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                        "model_path": model_path.as_posix(),
                    }
                )
            self._write_tsv(model_registry, ["family_category", "model_sha256", "model_path"], model_rows)
            panel_path = calibration.write_control_panel(
                root / "control_panel.tsv",
                calibration.build_control_panel(reference_manifest, seed_registry, model_registry),
            )

            outputs = calibration.build_calibration_command_manifest(panel_path, model_registry, root / "calibration")

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["family_category"] for row in rows], ["family_a", "family_b"])
            self.assertTrue(all("hmmsearch --noali --acc --seed 42 --cpu 1 --domtblout" in row["command"] for row in rows))
            self.assertTrue(all(row["command_status"] == "planned_not_run" for row in rows))
            first_fasta = Path(rows[0]["target_fasta_path"])
            self.assertIn(">family_a|cross_family_challenge|family_b|B1", first_fasta.read_text(encoding="ascii"))

    def test_calibration_commands_reject_control_panel_model_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_manifest, seed_registry, model_registry = self._build_fixture(root)
            panel_rows = calibration.build_control_panel(reference_manifest, seed_registry, model_registry)
            panel_rows[0]["model_sha256"] = "f" * 64
            panel_path = calibration.write_control_panel(root / "control_panel.tsv", panel_rows)
            model_path = root / "family_a.hmm"
            model_path.write_text("HMMER3/f family_a\n", encoding="ascii")
            self._write_tsv(
                model_registry,
                ["family_category", "model_sha256", "model_path"],
                [
                    {
                        "family_category": "family_a",
                        "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                        "model_path": model_path.as_posix(),
                    },
                    {
                        "family_category": "family_b",
                        "model_sha256": "b" * 64,
                        "model_path": model_path.as_posix(),
                    },
                ],
            )

            with self.assertRaisesRegex(ValueError, "model SHA256"):
                calibration.build_calibration_command_manifest(panel_path, model_registry, root / "calibration")


if __name__ == "__main__":
    unittest.main()
