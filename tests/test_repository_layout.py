"""Repository-level regression checks for the PHB-GTDB-GPT source tree."""

from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RepositoryLayoutTest(unittest.TestCase):
    def test_required_repository_files_exist(self) -> None:
        required = [
            "README.md",
            ".gitignore",
            ".gitattributes",
            "main.nf",
            "nextflow.config",
            "config/project.yaml",
            "config/paths.example.yaml",
            "docs/PROJECT_SCOPE.md",
            "scripts/validate_repository.py",
        ]
        missing = [path for path in required if not (ROOT / path).is_file()]
        self.assertEqual(missing, [])

    def test_large_input_directories_are_ignored(self) -> None:
        ignored = subprocess.run(
            ["git", "check-ignore", "00_raw_gtdb_r232/example.fna.gz"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(ignored.returncode, 0, ignored.stderr)

    def test_machine_specific_paths_are_ignored(self) -> None:
        ignored = subprocess.run(
            ["git", "check-ignore", "config/paths.yaml"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(ignored.returncode, 0, ignored.stderr)

    def test_generated_p06_outputs_are_ignored(self) -> None:
        paths = [
            "05_hmmer_scan/p06_hmmer_scan_manifest.tsv",
            "05_hmmer_scan/p06_hmmer_scan_summary.tsv",
            "05_hmmer_scan/p06_hmmer_candidates.tsv",
            "05_hmmer_scan/p06_hmmer_candidate_summary.tsv",
            "05_hmmer_scan/raw_domtblout/family/shard.domtblout",
            "05_hmmer_scan/hmmer_logs/family/shard.txt",
            "05_hmmer_scan/run_status/p06_hmmer_run_status.tsv",
            "05_hmmer_scan/overlong_protein_exclusions/family/chunk.tsv",
            "05_hmmer_scan_smoke_r8/raw_domtblout/family/shard.domtblout",
            "05_hmmer_scan_length_smoke_r8/overlong_protein_exclusions/family/chunk.tsv",
        ]
        for path in paths:
            with self.subTest(path=path):
                ignored = subprocess.run(
                    ["git", "check-ignore", path],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(ignored.returncode, 0, ignored.stderr)

    def test_repository_validator_succeeds(self) -> None:
        result = subprocess.run(
            [sys.executable, "scripts/validate_repository.py"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
