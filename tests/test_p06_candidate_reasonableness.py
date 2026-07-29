from __future__ import annotations

import csv
import tempfile
from pathlib import Path
import unittest

from scripts import p06_candidate_reasonableness as audit


class P06CandidateReasonablenessTest(unittest.TestCase):
    def _write_candidate_table(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "family_category",
                    "proteome_shard",
                    "target_id",
                    "full_sequence_score",
                    "tier",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "family_category": "extracellular_pha_depolymerase_core",
                    "proteome_shard": "chunk_000001",
                    "target_id": "protA",
                    "full_sequence_score": "190.0",
                    "tier": "High-confidence",
                }
            )
            writer.writerow(
                {
                    "family_category": "intracellular_phaZ_no_lipase_box",
                    "proteome_shard": "chunk_000001",
                    "target_id": "protA",
                    "full_sequence_score": "420.0",
                    "tier": "High-confidence",
                }
            )
            writer.writerow(
                {
                    "family_category": "intracellular_mcl_pha_dep",
                    "proteome_shard": "chunk_000002",
                    "target_id": "protB",
                    "full_sequence_score": "700.0",
                    "tier": "High-confidence",
                }
            )
            writer.writerow(
                {
                    "family_category": "archaeal_patatin_like_pha_dep",
                    "proteome_shard": "chunk_000003",
                    "target_id": "protC",
                    "full_sequence_score": "80.0",
                    "tier": "Review",
                }
            )

    def test_audit_writes_family_tier_unique_and_overlap_summaries(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table = root / "05_hmmer_scan" / "p06_hmmer_candidates.tsv"
            self._write_candidate_table(candidate_table)

            outputs = audit.audit_p06_candidates(
                candidate_table,
                root / "06_domain_annotation" / "p06_reasonableness",
                total_predicted_genes=1000,
                total_genomes=10,
            )

            self.assertTrue(outputs["summary"].is_file())
            self.assertTrue(outputs["family_tier_summary"].is_file())
            self.assertTrue(outputs["overlap_summary"].is_file())
            self.assertTrue(outputs["markdown_report"].is_file())

            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["family_tier_summary"].open("r", encoding="utf-8", newline="") as handle:
                family_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["overlap_summary"].open("r", encoding="utf-8", newline="") as handle:
                overlap_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertIn({"kind": "total", "name": "candidate_rows", "value": "4"}, summary_rows)
            self.assertIn({"kind": "total", "name": "unique_candidate_targets", "value": "3"}, summary_rows)
            self.assertIn({"kind": "tier", "name": "High-confidence_unique_targets", "value": "2"}, summary_rows)
            self.assertIn({"kind": "rate", "name": "High-confidence_rows_per_predicted_gene", "value": "0.003000"}, summary_rows)
            self.assertIn({"kind": "rate", "name": "High-confidence_unique_targets_per_genome", "value": "0.200000"}, summary_rows)
            self.assertIn(
                {
                    "family_category": "intracellular_phaZ_no_lipase_box",
                    "tier": "High-confidence",
                    "candidate_rows": "1",
                    "unique_targets": "1",
                    "proteome_shards": "1",
                    "max_full_sequence_score": "420.0",
                },
                family_rows,
            )
            self.assertEqual(overlap_rows[0]["target_id"], "protA")
            self.assertEqual(overlap_rows[0]["high_confidence_family_count"], "2")
            self.assertIn("not phenotype proof", outputs["markdown_report"].read_text(encoding="utf-8"))

    def test_audit_blocks_missing_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table = root / "bad.tsv"
            candidate_table.write_text("family_category\ttier\nx\tHigh-confidence\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "missing required columns"):
                audit.audit_p06_candidates(candidate_table, root / "audit")


if __name__ == "__main__":
    unittest.main()
