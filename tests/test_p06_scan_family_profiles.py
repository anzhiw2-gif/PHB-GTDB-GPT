from __future__ import annotations

import csv
import hashlib
import tempfile
from pathlib import Path
import unittest

from scripts import p06_scan_family_profiles as p06


class P06ScanFamilyProfilesTest(unittest.TestCase):
    def _write_inputs(self, root: Path) -> tuple[Path, Path, Path]:
        hmm_dir = root / "04_family_profiles" / "hmms"
        proteome_dir = root / "03_gtdb_proteomes" / "faa"
        hmm_dir.mkdir(parents=True, exist_ok=True)
        proteome_dir.mkdir(parents=True, exist_ok=True)

        (hmm_dir / "phaZ7_like.hmm").write_text("HMMER3/f dummy\n", encoding="utf-8")
        (hmm_dir / "intracellular_phaZ_no_lipase_box.hmm").write_text("HMMER3/f dummy\n", encoding="utf-8")
        (proteome_dir / "GCA_000001.faa").write_text(">prot1\nMMMMMMMMMM\n", encoding="utf-8")
        (proteome_dir / "GCA_000002.faa").write_text(">prot2\nMMMMMMMMMMMM\n", encoding="utf-8")
        registry = root / "04_family_profiles" / "manifests" / "p05_hmm_model_registry.tsv"
        registry.parent.mkdir(parents=True, exist_ok=True)
        self._write_model_registry(
            registry,
            [hmm_dir / "phaZ7_like.hmm", hmm_dir / "intracellular_phaZ_no_lipase_box.hmm"],
        )
        return hmm_dir, proteome_dir, registry

    def _write_model_registry(
        self,
        registry: Path,
        hmm_paths: list[Path],
        *,
        approved: bool = True,
        score_threshold: str = "120.0",
        hmm_coverage_threshold: str = "0.700000",
    ) -> None:
        with registry.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "family_category",
                    "approved_for_p06",
                    "scan_permission",
                    "model_path",
                    "model_sha256",
                    "model_specific_thresholds",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for hmm_path in hmm_paths:
                writer.writerow(
                    {
                        "family_category": hmm_path.stem,
                        "approved_for_p06": "yes" if approved else "no",
                        "scan_permission": "approved" if approved else "blocked",
                        "model_path": hmm_path.as_posix(),
                        "model_sha256": hashlib.sha256(hmm_path.read_bytes()).hexdigest(),
                        "model_specific_thresholds": (
                            f"full_score>={score_threshold};hmm_coverage>={hmm_coverage_threshold}" if approved else "none"
                        ),
                    }
                )

    def _write_domtblout(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "\n".join(
                [
                    "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias hmmfrom hmmto alifrom alito envfrom envto acc description",
                    "protA - 360 phaZ7_like - 400 1e-50 220.0 0.1 1 1 1e-60 1e-55 215.0 0.1 20 380 10 350 8 352 0.99 strong hit",
                    "protB - 380 phaZ7_like - 400 1e-08 45.0 0.4 1 1 1e-10 1e-09 43.0 0.4 40 200 30 170 28 172 0.88 review hit",
                    "protC - 220 phaZ7_like - 400 0.05 9.0 3.5 1 1 0.1 0.08 8.0 3.5 1 60 3 62 1 65 0.41 weak hit",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    def test_build_scan_manifest_writes_deterministic_command_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir, proteome_dir, registry = self._write_inputs(root)

            outputs = p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)

            self.assertTrue(outputs["manifest"].is_file())
            self.assertTrue(outputs["summary"].is_file())

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 4)
            self.assertEqual(
                [row["family_category"] for row in rows],
                ["intracellular_phaZ_no_lipase_box", "intracellular_phaZ_no_lipase_box", "phaZ7_like", "phaZ7_like"],
            )
            self.assertIn("hmmsearch --noali --acc --seed 42 --cpu 1 --domtblout", rows[0]["command"])
            self.assertEqual(len(rows[0]["model_sha256"]), 64)
            self.assertTrue(rows[0]["domtblout_path"].endswith(".domtblout"))
            self.assertTrue(rows[0]["overlong_exclusion_path"].endswith(".tsv"))
            self.assertIn("p06_stream_proteomes.py", rows[0]["command"])
            self.assertIn("set -o pipefail;", rows[0]["command"])
            self.assertEqual(rows[0]["command_status"], "planned_not_run")

    def test_build_scan_manifest_carries_calibrated_model_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir, proteome_dir, registry = self._write_inputs(root)
            self._write_model_registry(
                registry,
                [hmm_dir / "phaZ7_like.hmm", hmm_dir / "intracellular_phaZ_no_lipase_box.hmm"],
                score_threshold="300.0",
                hmm_coverage_threshold="0.800000",
            )

            outputs = p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(all(row["calibrated_full_score_threshold"] == "300.0" for row in rows))
            self.assertTrue(all(row["calibrated_hmm_coverage_threshold"] == "0.800000" for row in rows))

    def test_build_scan_manifest_accepts_recursive_gzipped_p03_proteomes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir = root / "04_family_profiles" / "hmms"
            proteome_dir = root / "03_gtdb_proteomes" / "faa"
            hmm_dir.mkdir(parents=True, exist_ok=True)
            nested_dir = proteome_dir / "GCF" / "001" / "234" / "567"
            nested_dir.mkdir(parents=True, exist_ok=True)

            (hmm_dir / "intracellular_mcl_pha_dep.hmm").write_text("HMMER3/f dummy\n", encoding="utf-8")
            (nested_dir / "GCF_001234567.1.faa.gz").write_text(">prot1\nMMMM\n", encoding="utf-8")

            registry = root / "04_family_profiles" / "manifests" / "p05_hmm_model_registry.tsv"
            registry.parent.mkdir(parents=True, exist_ok=True)
            self._write_model_registry(registry, [hmm_dir / "intracellular_mcl_pha_dep.hmm"])

            outputs = p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["proteome_shard"], "GCF_001234567.1")
            self.assertTrue(rows[0]["proteome_path"].endswith("GCF_001234567.1.faa.gz"))
            self.assertIn("GCF/001/234/567/GCF_001234567.1.faa.gz", rows[0]["command"])

    def test_build_scan_manifest_can_group_proteomes_into_streamed_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir, proteome_dir, registry = self._write_inputs(root)
            (proteome_dir / "GCA_000003.faa.gz").write_text(">prot3\nMMMMMMMM\n", encoding="utf-8")

            outputs = p06.build_scan_manifest(
                hmm_dir,
                proteome_dir,
                root / "05_hmmer_scan",
                model_registry=registry,
                proteomes_per_job=2,
            )

            with outputs["manifest"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[0]["proteome_shard"], "chunk_000001")
            self.assertEqual(rows[0]["proteome_count"], "2")
            self.assertIn("p06_stream_proteomes.py", rows[0]["command"])
            self.assertIn(" | hmmsearch --noali --acc --seed 42 --cpu 1 ", rows[0]["command"])
            self.assertTrue(rows[0]["command"].endswith(" -"))
            self.assertEqual(rows[1]["proteome_shard"], "chunk_000002")
            self.assertEqual(rows[1]["proteome_count"], "1")

    def test_parse_scan_manifest_classifies_candidates_from_domtblout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir, proteome_dir, registry = self._write_inputs(root)
            outputs = p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)
            manifest_path = outputs["manifest"]

            with manifest_path.open("r", encoding="utf-8", newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

            for row in manifest_rows:
                self._write_domtblout(Path(row["domtblout_path"]))

            parsed_outputs = p06.parse_scan_manifest(manifest_path, root / "05_hmmer_scan")

            self.assertTrue(parsed_outputs["candidates"].is_file())
            self.assertTrue(parsed_outputs["summary"].is_file())

            with parsed_outputs["candidates"].open("r", encoding="utf-8", newline="") as handle:
                candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
            with parsed_outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(candidate_rows), 12)
            tiers = {row["tier"] for row in candidate_rows}
            self.assertEqual(tiers, {"High-confidence", "Review", "Rejected"})

            strong = next(row for row in candidate_rows if row["target_id"] == "protA")
            review = next(row for row in candidate_rows if row["target_id"] == "protB")
            weak = next(row for row in candidate_rows if row["target_id"] == "protC")

            self.assertEqual(strong["tier"], "High-confidence")
            self.assertEqual(review["tier"], "Review")
            self.assertEqual(weak["tier"], "Rejected")
            self.assertIn("coverage", strong["tier_reason"])
            self.assertIn({"kind": "tier", "name": "High-confidence", "count": "4"}, summary_rows)
            self.assertIn({"kind": "tier", "name": "Review", "count": "4"}, summary_rows)
            self.assertIn({"kind": "tier", "name": "Rejected", "count": "4"}, summary_rows)

    def test_build_scan_manifest_rejects_unapproved_or_changed_models(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hmm_dir, proteome_dir, registry = self._write_inputs(root)
            self._write_model_registry(registry, [hmm_dir / "phaZ7_like.hmm"], approved=False)

            with self.assertRaisesRegex(ValueError, "no P06-approved HMMs"):
                p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)

            self._write_model_registry(registry, [hmm_dir / "phaZ7_like.hmm"])
            (hmm_dir / "phaZ7_like.hmm").write_text("changed model\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                p06.build_scan_manifest(hmm_dir, proteome_dir, root / "05_hmmer_scan", model_registry=registry)


if __name__ == "__main__":
    unittest.main()
