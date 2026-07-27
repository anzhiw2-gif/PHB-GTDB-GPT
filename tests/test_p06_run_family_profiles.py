from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path
import unittest

from scripts import p06_run_family_profiles as runner


class P06RunFamilyProfilesTest(unittest.TestCase):
    def _write_manifest(self, path: Path, domtblout_path: Path, command: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["family_category", "proteome_shard", "domtblout_path", "command"], delimiter="\t")
            writer.writeheader()
            writer.writerow(
                {
                    "family_category": "family_a",
                    "proteome_shard": "shard_1",
                    "domtblout_path": domtblout_path.as_posix(),
                    "command": command,
                }
            )

    def test_run_manifest_executes_missing_job_and_then_skips_completed_domtblout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            domtblout = root / "raw" / "family_a.domtblout"
            completer = root / "complete_domtblout.py"
            completer.write_text(
                "from pathlib import Path\n"
                "import sys\n"
                "path = Path(sys.argv[1])\n"
                "path.parent.mkdir(parents=True, exist_ok=True)\n"
                "path.write_text('# completed\\n', encoding='utf-8')\n",
                encoding="utf-8",
            )
            command = f'"{sys.executable}" "{completer}" "{domtblout}"'
            self._write_manifest(manifest, domtblout, command)

            first = runner.run_manifest(manifest, root / "status", workers=1)
            second = runner.run_manifest(manifest, root / "status", workers=1)

            self.assertEqual(first["completed"], 1)
            self.assertEqual(first["skipped_existing"], 0)
            self.assertEqual(second["completed"], 0)
            self.assertEqual(second["skipped_existing"], 1)
            self.assertTrue(domtblout.is_file())
            with (root / "status" / "p06_hmmer_run_status.tsv").open("r", encoding="utf-8", newline="") as handle:
                status_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(status_rows[0]["status"], "skipped_existing")

    def test_run_manifest_reexecutes_nonempty_output_without_a_completed_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            domtblout = root / "raw" / "family_a.domtblout"
            domtblout.parent.mkdir(parents=True, exist_ok=True)
            domtblout.write_text("# incomplete output\n", encoding="utf-8")
            completer = root / "complete_domtblout.py"
            completer.write_text(
                "from pathlib import Path\n"
                "import sys\n"
                "path = Path(sys.argv[1])\n"
                "path.write_text('# replaced complete output\\n', encoding='utf-8')\n",
                encoding="utf-8",
            )
            command = f'"{sys.executable}" "{completer}" "{domtblout}"'
            self._write_manifest(manifest, domtblout, command)

            result = runner.run_manifest(manifest, root / "status", workers=1)

            self.assertEqual(result["completed"], 1)
            self.assertEqual(result["skipped_existing"], 0)
            self.assertEqual(domtblout.read_text(encoding="utf-8"), "# replaced complete output\n")


if __name__ == "__main__":
    unittest.main()
