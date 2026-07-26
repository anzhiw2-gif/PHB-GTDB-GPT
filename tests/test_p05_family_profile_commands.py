from __future__ import annotations

import csv
import tempfile
from pathlib import Path
import unittest

from scripts import p05_family_profile_commands as commands
from scripts import p05_plan_family_profiles as p05


class P05FamilyProfileCommandsTest(unittest.TestCase):
    def _queue_row(self, **overrides: str) -> dict[str, str]:
        row = {field: "" for field in p05.BUILD_QUEUE_FIELDNAMES}
        row.update(
            {
                "family_category": "phaZ7_like",
                "taxonomic_domain": "Bacteria",
                "reference_library": "bacteria_high_confidence",
                "seed_row_count": "3",
                "qualifying_seed_row_count": "3",
                "independent_qualifying_accession_count": "3",
                "qualifying_seed_ids": "seed-a1;seed-a2;seed-a3",
                "qualifying_source_accessions": "A1;A2;A3",
                "seed_bundle_path": "seed_bundles/phaZ7_like.faa",
                "bundled_sequence_count": "3",
                "alignment_tool": "MAFFT",
                "alignment_mode": "L-INS-i",
                "hmm_build_tool": "hmmbuild",
                "calibration_search_tool": "hmmsearch",
                "calibration_control_panel": "close_non_target_hydrolases",
                "source_manifest_path": "reference_library.normalized.tsv",
                "source_plan_path": "p05_family_profile_plan.tsv",
                "notes": "synthetic eligible family",
            }
        )
        row.update(overrides)
        return row

    def _write_queue(self, path: Path, rows: list[dict[str, str]]) -> None:
        p05.write_tsv(path, rows, p05.BUILD_QUEUE_FIELDNAMES)

    def test_empty_queue_writes_header_only_manifest_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue_path = root / "p05_family_hmm_build_scaffold_queue.tsv"
            outdir = root / "manifests"
            self._write_queue(queue_path, [])

            outputs = commands.build_family_profile_commands(queue_path, outdir)

            self.assertEqual(
                outputs["manifest"].read_text(encoding="utf-8"),
                "\t".join(commands.COMMAND_FIELDNAMES) + "\n",
            )
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn({"kind": "total", "name": "scaffold_queue_rows", "count": "0"}, summary_rows)
            self.assertIn({"kind": "total", "name": "eligible_families", "count": "0"}, summary_rows)
            self.assertIn({"kind": "total", "name": "command_manifest_rows", "count": "0"}, summary_rows)

    def test_eligible_family_generates_mafft_and_hmmbuild_commands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue_path = root / "p05_family_hmm_build_scaffold_queue.tsv"
            outdir = root / "manifests"
            self._write_queue(queue_path, [self._queue_row(seed_bundle_path=r"seed_bundles\phaZ7_like.faa")])

            outputs = commands.build_family_profile_commands(queue_path, outdir)

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["family_category"], "phaZ7_like")
            self.assertEqual(row["seed_bundle_path"], "seed_bundles/phaZ7_like.faa")
            self.assertEqual(row["alignment_path"], (root / "alignments" / "phaZ7_like.aligned.faa").as_posix())
            self.assertEqual(row["hmm_path"], (root / "hmms" / "phaZ7_like.hmm").as_posix())
            self.assertIn("mafft --localpair --maxiterate 1000 --inputorder", row["alignment_command"])
            self.assertIn(" > ", row["alignment_command"])
            self.assertEqual(
                row["hmmbuild_command"],
                "hmmbuild --amino "
                f"'{(root / 'hmms' / 'phaZ7_like.hmm').as_posix()}' "
                f"'{(root / 'alignments' / 'phaZ7_like.aligned.faa').as_posix()}'",
            )
            self.assertEqual(row["command_status"], "planned_not_run")
            self.assertIn("not run", row["notes"])
            self.assertTrue((root / "alignments").is_dir())
            self.assertTrue((root / "hmms").is_dir())
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn({"kind": "total", "name": "eligible_families", "count": "1"}, summary_rows)

    def test_output_paths_and_commands_are_deterministic(self) -> None:
        rows = [
            self._queue_row(family_category="z_family", seed_bundle_path=r"bundles\z.faa"),
            self._queue_row(family_category="a_family", seed_bundle_path=r"bundles\a.faa"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue_path = root / "queue.tsv"
            self._write_queue(queue_path, rows)
            first = commands.build_family_profile_commands(queue_path, root / "manifests")
            first_content = first["manifest"].read_text(encoding="utf-8")

            second = commands.build_family_profile_commands(queue_path, root / "manifests")

            self.assertEqual(second["manifest"].read_text(encoding="utf-8"), first_content)
            with second["manifest"].open("r", encoding="utf-8", newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["family_category"] for row in output_rows], ["a_family", "z_family"])
            for row in output_rows:
                self.assertNotIn("\\", row["alignment_path"])
                self.assertNotIn("\\", row["hmm_path"])
                self.assertNotIn("\\", row["seed_bundle_path"])
                self.assertTrue(row["alignment_command"].endswith(f"'{row['alignment_path']}'"))
                self.assertTrue(row["hmmbuild_command"].endswith(f"'{row['alignment_path']}'"))

    def test_rejects_undersized_or_inconsistent_queue_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue_path = root / "p05_family_hmm_build_scaffold_queue.tsv"
            outdir = root / "manifests"
            cases = [
                ("at least 3", self._queue_row(independent_qualifying_accession_count="1", qualifying_source_accessions="A1")),
                (
                    "does not match independent_qualifying_accession_count",
                    self._queue_row(
                        independent_qualifying_accession_count="3",
                        qualifying_source_accessions="A1;A1;A2",
                    ),
                ),
            ]
            for expected_message, row in cases:
                self._write_queue(queue_path, [row])
                with self.assertRaisesRegex(ValueError, expected_message):
                    commands.build_family_profile_commands(queue_path, outdir)

    def test_rejects_bundle_count_inconsistency(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue_path = root / "p05_family_hmm_build_scaffold_queue.tsv"
            outdir = root / "manifests"
            self._write_queue(queue_path, [self._queue_row(bundled_sequence_count="2")])

            with self.assertRaisesRegex(ValueError, "bundled_sequence_count must match qualifying_seed_row_count"):
                commands.build_family_profile_commands(queue_path, outdir)


if __name__ == "__main__":
    unittest.main()
