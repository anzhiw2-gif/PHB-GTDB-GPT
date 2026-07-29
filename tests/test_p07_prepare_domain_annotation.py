from __future__ import annotations

import csv
import gzip
import tempfile
from pathlib import Path
import unittest
from unittest import mock

from scripts import p07_prepare_domain_annotation as p07


class P07PrepareDomainAnnotationTest(unittest.TestCase):
    def _write_fixture(self, root: Path) -> tuple[Path, Path]:
        proteome_dir = root / "03_gtdb_proteomes" / "faa"
        proteome_dir.mkdir(parents=True, exist_ok=True)
        proteome_path = proteome_dir / "GCA_000001.faa.gz"
        with gzip.open(proteome_path, "wt", encoding="utf-8") as handle:
            handle.write(">protA some P03 protein\n")
            handle.write("MMMMMMMMMM\n")
            handle.write(">protB\n")
            handle.write("ACDEFGHIKL\n")
            handle.write(">protC\n")
            handle.write("VVVVVVVVVVVV\n")

        scan_manifest = root / "05_hmmer_scan" / "p06_hmmer_scan_manifest.tsv"
        scan_manifest.parent.mkdir(parents=True, exist_ok=True)
        with scan_manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["family_category", "proteome_shard", "proteome_path"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "family_category": "intracellular_mcl_pha_dep",
                    "proteome_shard": "chunk_000001",
                    "proteome_path": proteome_path.as_posix(),
                }
            )
            writer.writerow(
                {
                    "family_category": "extracellular_pha_depolymerase_core",
                    "proteome_shard": "chunk_000001",
                    "proteome_path": proteome_path.as_posix(),
                }
            )

        candidate_table = root / "05_hmmer_scan" / "p06_hmmer_candidates.tsv"
        with candidate_table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "family_category",
                    "proteome_shard",
                    "target_id",
                    "target_length",
                    "full_sequence_score",
                    "tier",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "family_category": "intracellular_mcl_pha_dep",
                    "proteome_shard": "chunk_000001",
                    "target_id": "protA",
                    "target_length": "10",
                    "full_sequence_score": "220.0",
                    "tier": "High-confidence",
                }
            )
            writer.writerow(
                {
                    "family_category": "extracellular_pha_depolymerase_core",
                    "proteome_shard": "chunk_000001",
                    "target_id": "protA",
                    "target_length": "10",
                    "full_sequence_score": "180.0",
                    "tier": "High-confidence",
                }
            )
            writer.writerow(
                {
                    "family_category": "intracellular_mcl_pha_dep",
                    "proteome_shard": "chunk_000001",
                    "target_id": "protB",
                    "target_length": "10",
                    "full_sequence_score": "45.0",
                    "tier": "Review",
                }
            )
        return candidate_table, scan_manifest

    def test_prepare_p07_inputs_extracts_high_confidence_sequences_and_commands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table, scan_manifest = self._write_fixture(root)

            outputs = p07.prepare_p07_inputs(
                candidate_table,
                scan_manifest,
                root / "06_domain_annotation",
                sequences_per_shard=1,
                interproscan_cpu=2,
            )

            self.assertTrue(outputs["sequence_manifest"].is_file())
            self.assertTrue(outputs["command_manifest"].is_file())
            self.assertTrue(outputs["summary"].is_file())

            with outputs["sequence_manifest"].open("r", encoding="utf-8", newline="") as handle:
                sequence_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["command_manifest"].open("r", encoding="utf-8", newline="") as handle:
                command_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(sequence_rows), 1)
            self.assertEqual(sequence_rows[0]["target_id"], "protA")
            self.assertEqual(sequence_rows[0]["sequence_length"], "10")
            self.assertEqual(
                sequence_rows[0]["family_categories"],
                "extracellular_pha_depolymerase_core;intracellular_mcl_pha_dep",
            )
            self.assertEqual(sequence_rows[0]["max_full_sequence_score"], "220.0")
            self.assertIn("GTDB Release 11 R232", sequence_rows[0]["gtdb_release"])

            fasta_path = Path(sequence_rows[0]["fasta_shard"])
            self.assertTrue(fasta_path.is_file())
            fasta_text = fasta_path.read_text(encoding="utf-8")
            self.assertIn(">chunk_000001__protA", fasta_text)
            self.assertIn("MMMMMMMMMM", fasta_text)

            self.assertEqual([row["tool"] for row in command_rows], ["InterProScan", "SignalP6"])
            self.assertIn("-f TSV,JSON,GFF3", command_rows[0]["command"])
            self.assertIn("-cpu 2", command_rows[0]["command"])
            self.assertIn("--organism other", command_rows[1]["command"])
            self.assertTrue(all(row["command_status"] == "planned_not_run" for row in command_rows))
            self.assertIn(
                {"kind": "total", "name": "selected_unique_candidate_sequences", "count": "1"},
                summary_rows,
            )

    def test_prepare_p07_inputs_can_include_review_and_optional_phobius(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table, scan_manifest = self._write_fixture(root)

            outputs = p07.prepare_p07_inputs(
                candidate_table,
                scan_manifest,
                root / "06_domain_annotation",
                include_tiers=("High-confidence", "Review"),
                sequences_per_shard=2,
                phobius_exe="phobius.pl",
            )

            with outputs["sequence_manifest"].open("r", encoding="utf-8", newline="") as handle:
                sequence_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["command_manifest"].open("r", encoding="utf-8", newline="") as handle:
                command_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual([row["target_id"] for row in sequence_rows], ["protA", "protB"])
            self.assertEqual([row["tool"] for row in command_rows], ["InterProScan", "SignalP6", "Phobius"])
            self.assertIn("phobius.pl", command_rows[2]["command"])

    def test_prepare_p07_inputs_blocks_missing_selected_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table, scan_manifest = self._write_fixture(root)
            with candidate_table.open("a", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "family_category",
                        "proteome_shard",
                        "target_id",
                        "target_length",
                        "full_sequence_score",
                        "tier",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writerow(
                    {
                        "family_category": "archaeal_patatin_like_pha_dep",
                        "proteome_shard": "chunk_000001",
                        "target_id": "missing_protein",
                        "target_length": "10",
                        "full_sequence_score": "160.0",
                        "tier": "High-confidence",
                    }
                )

            with self.assertRaisesRegex(ValueError, "selected P06 candidate sequences were not found"):
                p07.prepare_p07_inputs(candidate_table, scan_manifest, root / "06_domain_annotation")

            missing_path = root / "06_domain_annotation" / "review" / "p07_missing_candidate_sequences.tsv"
            self.assertTrue(missing_path.is_file())

    def test_prepare_p07_inputs_blocks_length_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate_table, scan_manifest = self._write_fixture(root)
            text = candidate_table.read_text(encoding="utf-8")
            candidate_table.write_text(text.replace("protA\t10\t220.0", "protA\t99\t220.0"), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "target_length mismatch"):
                p07.prepare_p07_inputs(candidate_table, scan_manifest, root / "06_domain_annotation")

    def test_extract_records_uses_target_accession_to_avoid_unrelated_proteomes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wanted_path = root / "03_gtdb_proteomes" / "faa" / "GCA_000001.faa"
            unrelated_path = root / "03_gtdb_proteomes" / "faa" / "GCA_000002.faa"
            wanted_path.parent.mkdir(parents=True, exist_ok=True)
            wanted_path.write_text(">GCA_000001.1|contig|7\nMMMM\n", encoding="utf-8")
            unrelated_path.write_text(">GCA_000002.1|contig|9\nBAD\n", encoding="utf-8")
            groups = {
                ("chunk_000001", "GCA_000001.1|contig|7"): p07.CandidateGroup(
                    proteome_shard="chunk_000001",
                    target_id="GCA_000001.1|contig|7",
                )
            }
            groups[("chunk_000001", "GCA_000001.1|contig|7")].add_row(
                {
                    "family_category": "intracellular_phaZ_no_lipase_box",
                    "tier": "High-confidence",
                    "target_length": "4",
                    "full_sequence_score": "500.0",
                }
            )

            real_iter = p07.iter_fasta_records
            visited_paths: list[Path] = []

            def counting_iter(path: Path):
                visited_paths.append(path)
                yield from real_iter(path)

            with mock.patch.object(p07, "iter_fasta_records", side_effect=counting_iter):
                found, missing = p07.extract_records(
                    groups,
                    {"chunk_000001": [unrelated_path, wanted_path]},
                )

            self.assertEqual(missing, [])
            self.assertEqual(len(found), 1)
            self.assertEqual(found[0][1].source_path, wanted_path)
            self.assertEqual(visited_paths, [wanted_path])


if __name__ == "__main__":
    unittest.main()
