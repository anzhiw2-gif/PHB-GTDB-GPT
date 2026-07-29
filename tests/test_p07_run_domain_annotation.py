from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path
import unittest

from scripts import p07_run_domain_annotation as runner


class P07RunDomainAnnotationTest(unittest.TestCase):
    def _write_manifest(self, path: Path, rows: list[dict[str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "tool",
                    "fasta_shard",
                    "sequence_count",
                    "input_fasta",
                    "output_path",
                    "command",
                    "command_status",
                    "notes",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_preflight_manifest_reports_missing_executable_without_running_commands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / "should_not_exist.txt"
            manifest = root / "manifests" / "p07_domain_annotation_command_manifest.tsv"
            self._write_manifest(
                manifest,
                [
                    {
                        "tool": "InterProScan",
                        "fasta_shard": "p07_candidates_000001",
                        "sequence_count": "1",
                        "input_fasta": "input.faa",
                        "output_path": "interpro/p07_candidates_000001/interproscan",
                        "command": f"missing_interproscan_for_test --write {marker}",
                        "command_status": "planned_not_run",
                        "notes": "domain annotation only",
                    }
                ],
            )

            summary = runner.run_manifest(manifest, root / "run_status", workers=1, preflight_only=True)

            self.assertEqual(summary["missing_executable"], 1)
            self.assertFalse(marker.exists())
            with (root / "run_status" / "p07_domain_annotation_run_status.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                status_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(status_rows[0]["status"], "missing_executable")
            self.assertEqual(status_rows[0]["executable"], "missing_interproscan_for_test")

    def test_run_manifest_executes_job_and_then_skips_existing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifests" / "p07_domain_annotation_command_manifest.tsv"
            output_base = root / "interpro" / "p07_candidates_000001" / "interproscan"
            writer_script = root / "write_interpro_output.py"
            writer_script.write_text(
                "from pathlib import Path\n"
                "import sys\n"
                "base = Path(sys.argv[1])\n"
                "base.parent.mkdir(parents=True, exist_ok=True)\n"
                "base.with_suffix('.tsv').write_text('protein\\tPfam\\n', encoding='utf-8')\n",
                encoding="utf-8",
            )
            self._write_manifest(
                manifest,
                [
                    {
                        "tool": "InterProScan",
                        "fasta_shard": "p07_candidates_000001",
                        "sequence_count": "1",
                        "input_fasta": "input.faa",
                        "output_path": output_base.as_posix(),
                        "command": f'"{sys.executable}" "{writer_script}" "{output_base}"',
                        "command_status": "planned_not_run",
                        "notes": "domain annotation only",
                    }
                ],
            )

            first = runner.run_manifest(manifest, root / "run_status", workers=1)
            second = runner.run_manifest(manifest, root / "run_status", workers=1)

            self.assertEqual(first["completed"], 1)
            self.assertEqual(first["skipped_existing"], 0)
            self.assertEqual(second["completed"], 0)
            self.assertEqual(second["skipped_existing"], 1)
            self.assertTrue(output_base.with_suffix(".tsv").is_file())


if __name__ == "__main__":
    unittest.main()
