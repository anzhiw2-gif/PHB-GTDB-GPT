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

    def test_control_panel_sequence_hash_is_stable_across_fasta_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_manifest, seed_registry, model_registry = self._build_fixture(root)
            (root / "B1.faa").write_bytes(b">seed|B1\nMSTNPKPQRIT\n")
            first_rows = calibration.build_control_panel(reference_manifest, seed_registry, model_registry)
            first_hash = next(row["sequence_sha256"] for row in first_rows if row["source_accession"] == "B1")
            (root / "B1.faa").write_bytes(b">seed|B1\r\nMSTNPKPQRIT\r\n")

            second_rows = calibration.build_control_panel(reference_manifest, seed_registry, model_registry)
            second_hash = next(row["sequence_sha256"] for row in second_rows if row["source_accession"] == "B1")

            self.assertEqual(first_hash, second_hash)

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
            self.assertTrue(Path(rows[0]["domtblout_path"]).parent.is_dir())
            self.assertTrue(Path(rows[0]["main_output_path"]).parent.is_dir())

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

    def test_leave_one_out_commands_hold_out_each_seed_and_keep_three_training_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model_path = root / "family_a.hmm"
            model_path.write_text("HMMER3/f family_a\n", encoding="ascii")
            model_registry = root / "model_registry.tsv"
            self._write_tsv(
                model_registry,
                ["family_category", "model_sha256", "model_path"],
                [
                    {
                        "family_category": "family_a",
                        "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                        "model_path": model_path.as_posix(),
                    }
                ],
            )
            seed_rows: list[dict[str, str]] = []
            for accession in ("A1", "A2", "A3", "A4"):
                sequence_path = self._write_sequence(root, accession)
                seed_rows.append(
                    {
                        "family_category": "family_a",
                        "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
                        "seed_id": f"seed-{accession.lower()}",
                        "source_accession": accession,
                        "sequence_path": sequence_path.as_posix(),
                    }
                )
            seed_registry = root / "seed_registry.tsv"
            self._write_tsv(
                seed_registry,
                ["family_category", "model_sha256", "seed_id", "source_accession", "sequence_path"],
                seed_rows,
            )

            outputs = calibration.build_leave_one_out_command_manifest(seed_registry, model_registry, root / "calibration")

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["holdout_accession"] for row in rows], ["A1", "A2", "A3", "A4"])
            self.assertTrue(all(row["training_seed_count"] == "3" for row in rows))
            self.assertTrue(all("mafft --localpair --maxiterate 1000 --inputorder" in row["command"] for row in rows))
            self.assertTrue(all("hmmbuild --amino" in row["command"] for row in rows))
            first_bundle = Path(rows[0]["training_bundle_path"])
            self.assertNotIn("|A1\n", first_bundle.read_text(encoding="ascii"))
            self.assertIn(">positive|family_a|A1", Path(rows[0]["positive_fasta_path"]).read_text(encoding="ascii"))

    def test_parse_leave_one_out_results_reports_union_coverage_and_missing_positive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hit_domtblout = root / "hit.domtblout"
            hit_domtblout.write_text(
                "# domtblout fixture\n"
                "positive|family_a|A1 - 120 family_a - 100 1e-20 42.0 0.0 1 2 1e-20 1e-20 42.0 0.0 10 50 1 41 1 41 0.99 fixture\n"
                "positive|family_a|A1 - 120 family_a - 100 1e-20 42.0 0.0 2 2 1e-20 1e-20 10.0 0.0 45 80 42 77 42 77 0.99 fixture\n",
                encoding="ascii",
            )
            missing_domtblout = root / "missing.domtblout"
            missing_domtblout.write_text("# no reportable hit\n", encoding="ascii")
            manifest = root / "leave_one_out.tsv"
            self._write_tsv(
                manifest,
                list(calibration.LEAVE_ONE_OUT_COMMAND_FIELDNAMES),
                [
                    {
                        "family_category": "family_a",
                        "holdout_accession": "A1",
                        "domtblout_path": hit_domtblout.as_posix(),
                    },
                    {
                        "family_category": "family_b",
                        "holdout_accession": "B1",
                        "domtblout_path": missing_domtblout.as_posix(),
                    },
                ],
            )

            rows = calibration.parse_leave_one_out_results(manifest)

            self.assertEqual(rows[0]["positive_hit_status"], "recovered")
            self.assertEqual(rows[0]["best_full_score"], "42.0")
            self.assertEqual(rows[0]["hmm_coverage"], "0.710000")
            self.assertEqual(rows[0]["domain_count"], "2")
            self.assertEqual(rows[1]["positive_hit_status"], "missing")
            self.assertEqual(rows[1]["best_full_score"], "")
            self.assertEqual(rows[1]["hmm_coverage"], "")

    def test_derive_calibration_decisions_requires_positive_recovery_and_hard_challenge_separation(self) -> None:
        decisions = calibration.derive_calibration_decisions(
            [
                {"family_category": "family_a", "positive_hit_status": "recovered", "best_full_score": "42.0", "hmm_coverage": "0.700000"},
                {"family_category": "family_a", "positive_hit_status": "recovered", "best_full_score": "55.0", "hmm_coverage": "0.900000"},
                {"family_category": "family_b", "positive_hit_status": "missing", "best_full_score": "", "hmm_coverage": ""},
            ],
            [
                {"family_category": "family_a", "control_role": "cross_family_challenge", "hit_status": "hit", "best_full_score": "40.0", "hmm_coverage": "0.720000"},
                {"family_category": "family_a", "control_role": "cross_family_challenge", "hit_status": "hit", "best_full_score": "43.0", "hmm_coverage": "0.600000"},
                {"family_category": "family_b", "control_role": "cross_family_challenge", "hit_status": "no_hit", "best_full_score": "", "hmm_coverage": ""},
            ],
        )

        by_family = {row["family_category"]: row for row in decisions}
        self.assertEqual(by_family["family_a"]["proposed_score_threshold"], "42.0")
        self.assertEqual(by_family["family_a"]["proposed_hmm_coverage_threshold"], "0.700000")
        self.assertEqual(by_family["family_a"]["hard_challenges_passing_proposed_rule"], "0")
        self.assertEqual(by_family["family_a"]["recommendation"], "eligible_for_human_review")
        self.assertEqual(by_family["family_b"]["positive_recovery_missing"], "1")
        self.assertEqual(by_family["family_b"]["recommendation"], "blocked_positive_recovery_failed")

    def test_derive_calibration_decisions_rejects_family_without_hard_challenge(self) -> None:
        with self.assertRaisesRegex(ValueError, "hard cross-family challenge"):
            calibration.derive_calibration_decisions(
                [
                    {
                        "family_category": "family_a",
                        "positive_hit_status": "recovered",
                        "best_full_score": "42.0",
                        "hmm_coverage": "0.700000",
                    }
                ],
                [],
            )

    def test_parse_control_smoke_results_retains_unhit_controls(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            domtblout = root / "family_a.controls.domtblout"
            hard_control_id = "family_a|cross_family_challenge|family_b|B1"
            boundary_control_id = "family_a|boundary_observation|family_a|A2"
            domtblout.write_text(
                "# domtblout fixture\n"
                f"{hard_control_id} - 120 family_a - 100 1e-20 42.0 0.0 1 1 1e-20 1e-20 42.0 0.0 10 60 1 51 1 51 0.99 fixture\n",
                encoding="ascii",
            )
            panel = root / "panel.tsv"
            self._write_tsv(
                panel,
                list(calibration.CONTROL_PANEL_FIELDNAMES),
                [
                    {"family_category": "family_a", "control_id": hard_control_id, "control_role": "cross_family_challenge"},
                    {"family_category": "family_a", "control_id": boundary_control_id, "control_role": "boundary_observation"},
                ],
            )
            manifest = root / "control_commands.tsv"
            self._write_tsv(
                manifest,
                list(calibration.CALIBRATION_COMMAND_FIELDNAMES),
                [{"family_category": "family_a", "domtblout_path": domtblout.as_posix()}],
            )

            rows = calibration.parse_control_smoke_results(panel, manifest)

            by_id = {row["control_id"]: row for row in rows}
            self.assertEqual(by_id[hard_control_id]["hit_status"], "hit")
            self.assertEqual(by_id[hard_control_id]["best_full_score"], "42.0")
            self.assertEqual(by_id[hard_control_id]["hmm_coverage"], "0.510000")
            self.assertEqual(by_id[boundary_control_id]["hit_status"], "no_hit")

    def test_write_calibration_result_tables_writes_detailed_and_decision_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            outputs = calibration.write_calibration_result_tables(
                root,
                [
                    {
                        "family_category": "family_a",
                        "holdout_accession": "A1",
                        "domtblout_path": "ignored/raw.domtblout",
                        "positive_hit_status": "recovered",
                        "best_full_score": "42.0",
                        "hmm_coverage": "0.700000",
                        "domain_count": "1",
                    }
                ],
                [
                    {
                        "family_category": "family_a",
                        "control_id": "family_a|cross_family_challenge|family_b|B1",
                        "control_role": "cross_family_challenge",
                        "hit_status": "no_hit",
                        "best_full_score": "",
                        "hmm_coverage": "",
                        "domain_count": "0",
                    }
                ],
            )

            self.assertTrue(outputs["leave_one_out"].is_file())
            self.assertTrue(outputs["control_smoke"].is_file())
            with outputs["decisions"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["recommendation"], "eligible_for_human_review")


if __name__ == "__main__":
    unittest.main()
