from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

from scripts import p08_run_phylogeny as runner


MANIFEST_FIELDS = (
    "family_category", "command_status", "input_fasta_path", "input_sha256",
    "candidate_input_record_count", "total_input_record_count", "route", "alignment_fasta_path",
    "representative_input_fasta_path", "fasttree_tree_path", "iqtree_prefix", "representative_plan",
    "mafft_template", "fasttree_template", "iqtree2_template", "iqtree2_annotation",
    "rooting_policy", "evidence_boundary",
)


class P08RunPhylogenyTest(unittest.TestCase):
    def _write_fasta(self, path: Path, text: str = ">fixture\nMPEPTIDE\n") -> str:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def _row(
        self,
        root: Path,
        *,
        family: str = "family_a",
        route: str = "mafft_linsi_then_review",
        input_path: Path | None = None,
        input_sha256: str | None = None,
        output_path: Path | None = None,
        command: str = "",
    ) -> dict[str, str]:
        input_path = input_path or root / f"{family}.faa"
        input_sha256 = input_sha256 or self._write_fasta(input_path)
        output_path = output_path or root / "outputs" / f"{family}.aln.faa"
        return {
            "family_category": family,
            "command_status": "planned_not_run",
            "input_fasta_path": str(input_path),
            "input_sha256": input_sha256,
            "candidate_input_record_count": "1",
            "total_input_record_count": "1",
            "route": route,
            "alignment_fasta_path": str(output_path),
            "representative_input_fasta_path": str(input_path),
            "fasttree_tree_path": str(root / "outputs" / f"{family}.nwk"),
            "iqtree_prefix": str(root / "outputs" / family),
            "representative_plan": "",
            "mafft_template": command,
            "fasttree_template": command if "fasttree" in route else "",
            "iqtree2_template": "iqtree2 -s {alignment_fasta} --prefix {iqtree_prefix}",
            "iqtree2_annotation": "requires_independent_subset_and_outgroup_approval",
            "rooting_policy": "explicit_accessioned_outgroup_required; otherwise midpoint_display_only",
            "evidence_boundary": "sequence_and_annotation_evidence_only_not_phenotype_proof",
        }

    def _write_manifest(self, root: Path, rows: list[dict[str, str]]) -> Path:
        manifest = root / "manifests" / "p08_phylogeny_command_manifest.tsv"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return manifest

    def _read_status(self, status_dir: Path) -> list[dict[str, str]]:
        with (status_dir / runner.STATUS_FILENAME).open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def _writer_script(self, root: Path) -> Path:
        script = root / "writer.py"
        script.write_text(
            "from pathlib import Path\n"
            "import sys\n"
            "mode, output = sys.argv[1:]\n"
            "if mode == 'fail':\n"
            "    raise SystemExit(7)\n"
            "if mode == 'write':\n"
            "    path = Path(output)\n"
            "    path.parent.mkdir(parents=True, exist_ok=True)\n"
            "    path.write_text('tree planning fixture\\n', encoding='utf-8')\n",
            encoding="utf-8",
        )
        return script

    def test_preflight_missing_executable_does_not_execute_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / "must_not_exist.txt"
            manifest = self._write_manifest(root, [self._row(root, command=f"p08_missing_executable --marker {marker}")])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            self.assertEqual(summary["missing_executable"], 1)
            self.assertFalse(marker.exists())
            self.assertEqual(self._read_status(root / "status")[0]["status"], "missing_executable")

    def test_preflight_reports_missing_input_before_other_checks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "absent.faa"
            manifest = self._write_manifest(root, [self._row(root, input_path=missing, input_sha256="a" * 64, command="p08_missing")])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            self.assertEqual(summary["missing_input"], 1)

    def test_missing_input_precedes_malformed_command_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "absent.faa"
            manifest = self._write_manifest(root, [self._row(root, input_path=missing, input_sha256="a" * 64, command="'unterminated")])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            self.assertEqual(summary["missing_input"], 1)

    def test_preflight_reports_checksum_mismatch_before_command_checks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_path = root / "input.faa"
            self._write_fasta(input_path)
            manifest = self._write_manifest(root, [self._row(root, input_path=input_path, input_sha256="0" * 64, command="p08_missing")])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            self.assertEqual(summary["checksum_mismatch"], 1)

    def test_malformed_selected_command_is_a_manifest_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = self._write_manifest(root, [self._row(root, command="'unterminated")])

            with self.assertRaisesRegex(ValueError, "could not parse"):
                runner.run_manifest(manifest, root / "status", preflight_only=True)

    def test_preflight_safe_executable_does_not_create_planned_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = self._writer_script(root)
            output = root / "outputs" / "planned.aln.faa"
            command = f'"{sys.executable}" "{writer}" write "{output}"'
            manifest = self._write_manifest(root, [self._row(root, output_path=output, command=command)])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            rows = self._read_status(root / "status")
            self.assertEqual(summary["preflight_ok"], 1)
            self.assertFalse(output.exists())
            self.assertEqual(rows[0]["status"], "preflight_ok")
            self.assertIn("not biological negative evidence", rows[0]["notes"])

    def test_large_fasttree_preflight_verifies_full_family_input_before_representative_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_fasta = root / "family_full.faa"
            family_sha256 = self._write_fasta(family_fasta, ">full\nMPEPTIDE\n")
            representative_fasta = root / "representatives" / "family_representatives.faa"
            output = root / "trees" / "family.fasttree.nwk"
            row = self._row(
                root,
                family="large_family",
                route="deterministic_representative_plan_then_fasttree_exploratory",
                input_path=family_fasta,
                input_sha256=family_sha256,
                output_path=output,
                command="FastTree -lg {representative_alignment_fasta} > {fasttree_tree}",
            )
            row["representative_input_fasta_path"] = str(representative_fasta)
            manifest = self._write_manifest(root, [row])

            with patch("scripts.p08_run_phylogeny._executable_available", return_value=True):
                summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            status = self._read_status(root / "status")[0]
            self.assertEqual(summary["preflight_ok"], 1)
            self.assertFalse(representative_fasta.exists())
            self.assertEqual(status["input_fasta_path"], str(family_fasta))
            self.assertEqual(status["input_sha256"], family_sha256)
            self.assertEqual(status["representative_input_fasta_path"], str(representative_fasta))
            self.assertEqual(status["representative_input_sha256"], "")
            self.assertIn("independent representative SHA-256", status["representative_input_contract"])
            self.assertIn(str(representative_fasta), status["selected_command"])
            self.assertNotIn(str(family_fasta), status["selected_command"])

    def test_manifest_requires_planned_not_run_command_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = self._row(root, command="p08_missing")
            row["command_status"] = "completed"
            manifest = self._write_manifest(root, [row])

            with self.assertRaisesRegex(ValueError, "command_status=planned_not_run"):
                runner.run_manifest(manifest, root / "status", preflight_only=True)

    def test_existing_nonempty_route_output_is_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = self._writer_script(root)
            output = root / "outputs" / "already.aln.faa"
            output.parent.mkdir(parents=True)
            output.write_text("existing output\n", encoding="utf-8")
            command = f'"{sys.executable}" "{writer}" write "{output}"'
            manifest = self._write_manifest(root, [self._row(root, output_path=output, command=command)])

            summary = runner.run_manifest(manifest, root / "status", preflight_only=True)

            self.assertEqual(summary["skipped_existing"], 1)

    def test_status_file_is_atomic_sorted_and_replaces_interrupted_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = self._writer_script(root)
            command = f'"{sys.executable}" "{writer}" noop ignored'
            rows = [self._row(root, family="z_family", command=command), self._row(root, family="a_family", command=command)]
            manifest = self._write_manifest(root, rows)
            status_dir = root / "status"
            status_dir.mkdir()
            (status_dir / runner.STATUS_FILENAME).write_text("interrupted", encoding="utf-8")

            summary = runner.run_manifest(manifest, status_dir, workers=2, preflight_only=True)

            status_rows = self._read_status(status_dir)
            self.assertEqual(summary["preflight_ok"], 2)
            self.assertEqual([row["family_category"] for row in status_rows], ["a_family", "z_family"])
            self.assertEqual(list(status_dir.glob(f".{runner.STATUS_FILENAME}.*.tmp")), [])

    def test_safe_nonpreflight_writer_completes_then_resumes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = self._writer_script(root)
            output = root / "outputs" / "completed.aln.faa"
            command = f'"{sys.executable}" "{writer}" write "{output}"'
            manifest = self._write_manifest(root, [self._row(root, output_path=output, command=command)])

            first = runner.run_manifest(manifest, root / "status", preflight_only=False, allow_test_execution=True)
            second = runner.run_manifest(manifest, root / "status", preflight_only=False, allow_test_execution=True)

            self.assertEqual(first["completed"], 1)
            self.assertEqual(second["skipped_existing"], 1)
            self.assertTrue(output.is_file())

    def test_nonpreflight_writer_failure_and_missing_output_are_distinguished(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            writer = self._writer_script(root)
            failed_output = root / "outputs" / "failed.aln.faa"
            missing_output = root / "outputs" / "missing.aln.faa"
            failed = self._row(root, family="failed", output_path=failed_output, command=f'"{sys.executable}" "{writer}" fail "{failed_output}"')
            no_output = self._row(root, family="no_output", output_path=missing_output, command=f'"{sys.executable}" "{writer}" noop "{missing_output}"')
            manifest = self._write_manifest(root, [failed, no_output])

            summary = runner.run_manifest(manifest, root / "status", workers=2, preflight_only=False, allow_test_execution=True)

            self.assertEqual(summary["failed_exit_code"], 1)
            self.assertEqual(summary["failed_missing_output"], 1)

    def test_nonpreflight_rejects_real_phylogeny_executables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for executable in ("mafft", "iqtree2", "FastTree", "mafft.exe", "FastTree.exe", "iqtree2.exe", str(root / "bin" / "FastTree.exe")):
                with self.subTest(executable=executable):
                    manifest = self._write_manifest(root, [self._row(root, command=f"{executable} --version")])
                    with patch("scripts.p08_run_phylogeny._executable_available", return_value=True):
                        with patch("scripts.p08_run_phylogeny.subprocess.run") as run:
                            with self.assertRaisesRegex(ValueError, "not authorized"):
                                runner.run_manifest(manifest, root / "status", preflight_only=False, allow_test_execution=True)
                    run.assert_not_called()

    def test_nonpreflight_rejects_non_fixture_python_before_subprocess(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = self._write_manifest(root, [self._row(root, command=f'"{sys.executable}" -c "print(1)"')])

            with patch("scripts.p08_run_phylogeny.subprocess.run") as run:
                with self.assertRaisesRegex(ValueError, "test-fixture"):
                    runner.run_manifest(manifest, root / "status", preflight_only=False, allow_test_execution=True)

            run.assert_not_called()

    def test_cli_preflight_only_succeeds_without_executing_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / "cli_marker.txt"
            manifest = self._write_manifest(root, [self._row(root, command=f"p08_cli_missing --marker {marker}")])

            completed = subprocess.run(
                [sys.executable, "scripts/p08_run_phylogeny.py", "--manifest", str(manifest), "--status-dir", str(root / "status"), "--workers", "2", "--preflight-only"],
                cwd=Path(__file__).resolve().parents[1],
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("missing_executable: 1", completed.stdout)
            self.assertFalse(marker.exists())


if __name__ == "__main__":
    unittest.main()
