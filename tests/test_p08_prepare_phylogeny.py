"""Tests for P08 candidate, taxonomy, and reference input preparation."""

from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.p08_prepare_phylogeny import (
    DEFAULT_INCLUDE_TIERS,
    REQUIRED_P07_TOOLS,
    prepare_p08_inputs,
    route_family_size,
)


P06_FIELDS = (
    "family_category",
    "proteome_shard",
    "target_id",
    "target_accession",
    "target_length",
    "full_sequence_score",
    "hmm_coverage",
    "tier",
)
P07_SEQUENCE_FIELDS = (
    "p07_sequence_id",
    "proteome_shard",
    "target_id",
    "source_proteome_path",
    "target_length_from_p06",
    "sequence_length",
    "family_categories",
    "fasta_shard",
)


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PrepareP08InputsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.outdir = self.root / "p08"
        self.p06 = self.root / "p06.tsv"
        self.p07_sequences = self.root / "p07_sequences.tsv"
        self.p07_status = self.root / "p07_status.tsv"
        self.registry = self.root / "registry.tsv"
        self.seeds = self.root / "seeds.tsv"
        self.controls = self.root / "controls.tsv"
        self.taxonomy = self.root / "taxonomy.tsv"
        self.seed_bacterial = self._sequence_file("seed_bacterial.faa", ">seed-b\nMPEPTIDE\n")
        self.control_bacterial = self._sequence_file("control_bacterial.faa", ">control-b\nAAAA\n")
        self.seed_archaeal = self._sequence_file("seed_archaeal.faa", ">seed-a\nMKKK\n")
        self.control_archaeal = self._sequence_file("control_archaeal.faa", ">control-a\nVVVV\n")
        self._write_complete_fixture()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _sequence_file(self, name: str, content: str) -> Path:
        path = self.root / name
        path.write_text(content, encoding="utf-8")
        return path

    @staticmethod
    def _sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def _write_complete_fixture(self) -> None:
        write_tsv(
            self.p06,
            P06_FIELDS,
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "proteome_shard": "part_b", "target_id": "arc1", "target_accession": "arc-accession", "target_length": "4", "full_sequence_score": "200.0", "hmm_coverage": "0.8", "tier": "High-confidence"},
                {"family_category": "extracellular_pha_depolymerase_core", "proteome_shard": "part_a", "target_id": "bac1", "target_accession": "bac-accession", "target_length": "8", "full_sequence_score": "180.0", "hmm_coverage": "0.7", "tier": "High-confidence"},
                {"family_category": "extracellular_pha_depolymerase_core", "proteome_shard": "part_a", "target_id": "review1", "target_accession": "review-accession", "target_length": "6", "full_sequence_score": "160.0", "hmm_coverage": "0.5", "tier": "Review"},
            ],
        )
        write_tsv(
            self.p07_sequences,
            P07_SEQUENCE_FIELDS,
            [
                {"p07_sequence_id": "p07-arc1", "proteome_shard": "part_b", "target_id": "arc1", "source_proteome_path": "/machine/RS_GCF_000002.faa.gz", "target_length_from_p06": "4", "sequence_length": "4", "family_categories": "archaeal_patatin_like_pha_dep", "fasta_shard": "/machine/p07_arc.fa"},
                {"p07_sequence_id": "p07-bac1", "proteome_shard": "part_a", "target_id": "bac1", "source_proteome_path": "/machine/GB_GCF_000001.faa.gz", "target_length_from_p06": "8", "sequence_length": "8", "family_categories": "extracellular_pha_depolymerase_core", "fasta_shard": "/machine/p07_bac.fa"},
                {"p07_sequence_id": "p07-review1", "proteome_shard": "part_a", "target_id": "review1", "source_proteome_path": "/machine/GB_GCF_000001.faa.gz", "target_length_from_p06": "6", "sequence_length": "6", "family_categories": "extracellular_pha_depolymerase_core", "fasta_shard": "/machine/p07_bac.fa"},
            ],
        )
        write_tsv(
            self.p07_status,
            ("tool", "fasta_shard", "status"),
            [
                {"tool": tool, "fasta_shard": shard, "status": "completed"}
                for tool in REQUIRED_P07_TOOLS
                for shard in ("/machine/p07_arc.fa", "/machine/p07_bac.fa")
            ],
        )
        write_tsv(
            self.registry,
            ("family_category", "approved_for_p06", "scan_permission", "model_sha256"),
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "approved_for_p06": "yes", "scan_permission": "approved", "model_sha256": "a" * 64},
                {"family_category": "extracellular_pha_depolymerase_core", "approved_for_p06": "yes", "scan_permission": "approved", "model_sha256": "b" * 64},
            ],
        )
        write_tsv(
            self.seeds,
            ("family_category", "seed_id", "source_accession", "sequence_path", "sequence_sha256", "evidence", "notes"),
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "seed_id": "seed-arc", "source_accession": "ARC1", "sequence_path": str(self.seed_archaeal), "sequence_sha256": self._sha256(self.seed_archaeal), "evidence": "E3", "notes": "archaeal seed"},
                {"family_category": "extracellular_pha_depolymerase_core", "seed_id": "seed-bac", "source_accession": "BAC1", "sequence_path": str(self.seed_bacterial), "sequence_sha256": self._sha256(self.seed_bacterial), "evidence": "E1", "notes": "bacterial seed"},
            ],
        )
        write_tsv(
            self.controls,
            ("family_category", "control_id", "control_role", "sequence_path", "sequence_sha256", "evidence", "notes"),
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "control_id": "control-arc", "control_role": "hard_negative", "sequence_path": str(self.control_archaeal), "sequence_sha256": self._sha256(self.control_archaeal), "evidence": "control", "notes": "archaeal control"},
                {"family_category": "extracellular_pha_depolymerase_core", "control_id": "control-bac", "control_role": "hard_negative", "sequence_path": str(self.control_bacterial), "sequence_sha256": self._sha256(self.control_bacterial), "evidence": "control", "notes": "bacterial control"},
            ],
        )
        self.taxonomy.write_text("GCF_000001\td__Bacteria;p__Test\nGCF_000002\td__Archaea;p__Test\n", encoding="utf-8")

    def _prepare(self, **kwargs: object) -> dict[str, Path]:
        outdir = kwargs.pop("outdir", self.outdir)
        return prepare_p08_inputs(
            p06_candidate_table=self.p06,
            p07_sequence_table=self.p07_sequences,
            p07_status_table=self.p07_status,
            p05_model_registry=self.registry,
            p05_seed_table=self.seeds,
            p05_control_table=self.controls,
            taxonomy_paths=[self.taxonomy],
            outdir=outdir,
            **kwargs,
        )

    def test_high_confidence_candidates_have_sorted_taxonomy_annotations_and_verified_references(self) -> None:
        outputs = self._prepare()
        candidates = read_tsv(outputs["candidate_manifest"])
        taxonomy = read_tsv(outputs["taxonomy_join"])
        references = read_tsv(outputs["family_reference_manifest"])
        self.assertEqual(DEFAULT_INCLUDE_TIERS, ("High-confidence",))
        self.assertEqual([row["family_category"] for row in candidates], ["archaeal_patatin_like_pha_dep", "extracellular_pha_depolymerase_core"])
        self.assertEqual([row["taxonomy_lineage"] for row in taxonomy], ["d__Archaea;p__Test", "d__Bacteria;p__Test"])
        self.assertTrue(all(row["p07_annotation_status"] == "completed" for row in candidates))
        self.assertTrue(all(row["verified_sha256"] == row["sequence_sha256"] for row in references))
        self.assertTrue(all(row["evidence_boundary"] == "sequence_and_annotation_evidence_only_not_phenotype_proof" for row in candidates))

    def test_review_candidate_is_excluded_by_default_and_included_only_when_requested(self) -> None:
        self._prepare()
        default_rows = read_tsv(self.outdir / "manifests" / "p08_candidate_manifest.tsv")
        self.assertNotIn("review1", [row["target_id"] for row in default_rows])
        requested_outdir = self.root / "with-review"
        self._prepare(outdir=requested_outdir, include_tiers=("High-confidence", "Review"))
        requested_rows = read_tsv(requested_outdir / "manifests" / "p08_candidate_manifest.tsv")
        self.assertIn("review1", [row["target_id"] for row in requested_rows])

    def test_length_mismatch_blocks_and_writes_review(self) -> None:
        rows = read_tsv(self.p07_sequences)
        rows[0]["sequence_length"] = "5"
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, rows)
        with self.assertRaisesRegex(ValueError, "P06/P07 sequence length mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any("P06/P07 sequence length mismatch" in row["reason"] for row in blocks))

    def test_missing_taxonomy_blocks_and_writes_review(self) -> None:
        self.taxonomy.write_text("GCF_000001\td__Bacteria;p__Test\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "taxonomy"):
            self._prepare()
        self.assertTrue((self.outdir / "review" / "p08_blocked_records.tsv").exists())

    def test_missing_approved_registry_row_blocks_processing(self) -> None:
        rows = read_tsv(self.registry)
        rows[1]["approved_for_p06"] = "no"
        write_tsv(self.registry, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "approved P05 model"):
            self._prepare()

    def test_missing_p07_tool_status_blocks_with_candidate_and_status_details(self) -> None:
        rows = [row for row in read_tsv(self.p07_status) if not (row["tool"] == "InterProScan" and row["fasta_shard"] == "/machine/p07_arc.fa")]
        write_tsv(self.p07_status, ("tool", "fasta_shard", "status"), rows)
        with self.assertRaisesRegex(ValueError, "P07 annotation status requirement failed"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        block = next(row for row in blocks if row["target_id"] == "arc1")
        self.assertEqual(block["family_category"], "archaeal_patatin_like_pha_dep")
        self.assertEqual(block["proteome_shard"], "part_b")
        self.assertIn("InterProScan=missing", block["notes"])

    def test_failed_exit_code_p07_tool_status_blocks(self) -> None:
        rows = read_tsv(self.p07_status)
        rows[0]["status"] = "failed_exit_code"
        write_tsv(self.p07_status, ("tool", "fasta_shard", "status"), rows)
        with self.assertRaisesRegex(ValueError, "P07 annotation status requirement failed"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any("failed_exit_code" in row["notes"] for row in blocks))

    def test_planned_not_run_p07_tool_status_blocks(self) -> None:
        rows = read_tsv(self.p07_status)
        rows[0]["status"] = "planned_not_run"
        write_tsv(self.p07_status, ("tool", "fasta_shard", "status"), rows)
        with self.assertRaisesRegex(ValueError, "P07 annotation status requirement failed"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any("planned_not_run" in row["notes"] for row in blocks))

    def test_unknown_review_family_blocks_before_tier_selection(self) -> None:
        rows = read_tsv(self.p06)
        rows[2]["family_category"] = "unknown_review_family"
        write_tsv(self.p06, P06_FIELDS, rows)
        with self.assertRaisesRegex(ValueError, "unknown P06 family"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertEqual(blocks[0]["target_id"], "review1")

    def test_blocked_review_family_blocks_before_tier_selection(self) -> None:
        rows = read_tsv(self.p06)
        rows[2]["family_category"] = "intracellular_mcl_pha_dep"
        write_tsv(self.p06, P06_FIELDS, rows)
        registry_rows = read_tsv(self.registry)
        registry_rows.append({"family_category": "intracellular_mcl_pha_dep", "approved_for_p06": "no", "scan_permission": "blocked", "model_sha256": "c" * 64})
        write_tsv(self.registry, tuple(registry_rows[0].keys()), registry_rows)
        with self.assertRaisesRegex(ValueError, "missing approved P05 model"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertEqual(blocks[0]["target_id"], "review1")

    def test_seed_or_control_checksum_mismatch_blocks_and_writes_review(self) -> None:
        rows = read_tsv(self.seeds)
        rows[0]["sequence_sha256"] = "0" * 64
        write_tsv(self.seeds, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self._prepare()
        self.assertTrue((self.outdir / "review" / "p08_blocked_records.tsv").exists())

    def test_route_family_size_boundaries(self) -> None:
        self.assertEqual(route_family_size(199), "mafft_linsi_then_review")
        self.assertEqual(route_family_size(200), "mafft_auto_then_review")
        self.assertEqual(route_family_size(2001), "deterministic_representative_plan_then_fasttree_exploratory")

    def test_prepare_p08_inputs_never_calls_subprocess_run(self) -> None:
        with patch("scripts.p08_prepare_phylogeny.subprocess.run") as run:
            self._prepare()
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
