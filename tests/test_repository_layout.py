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

    def test_generated_p07_outputs_are_ignored(self) -> None:
        paths = [
            "06_domain_annotation/p06_reasonableness/p06_candidate_reasonableness_summary.tsv",
            "06_domain_annotation/p06_reasonableness/p06_family_tier_reasonableness.tsv",
            "06_domain_annotation/p06_reasonableness/p06_high_confidence_overlap_targets.tsv",
            "06_domain_annotation/p06_reasonableness/P06_REASONABLENESS_AUDIT.md",
            "06_domain_annotation/input/fasta_shards/p07_candidates_000001.faa",
            "06_domain_annotation/interpro/p07_candidates_000001/interproscan.tsv",
            "06_domain_annotation/localization/signalp6/p07_candidates_000001/prediction_results.txt",
            "06_domain_annotation/review/p07_missing_candidate_sequences.tsv",
            "06_domain_annotation/run_status/p07_domain_annotation_run_status.tsv",
            "06_domain_annotation/manifests/p07_candidate_sequence_manifest.tsv",
            "06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv",
            "06_domain_annotation/manifests/p07_domain_annotation_summary.tsv",
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

    def test_generated_p08_outputs_are_ignored(self) -> None:
        paths = [
            "07_phylogeny/manifests/p08_candidate_manifest.tsv",
            "07_phylogeny/manifests/p08_phylogeny_command_manifest.tsv",
            "07_phylogeny/review/p08_blocked_records.tsv",
            "07_phylogeny/gtdb_mapping/p08_taxonomy_join.tsv",
            "07_phylogeny/run_status/p08_phylogeny_run_status.tsv",
            "07_phylogeny/run_logs/mafft__family.stderr.log",
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
