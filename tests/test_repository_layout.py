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
