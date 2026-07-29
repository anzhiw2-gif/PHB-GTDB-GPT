"""Tests for P08 candidate, taxonomy, and reference input preparation."""

from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.p08_prepare_phylogeny import (
    DEFAULT_INCLUDE_TIERS,
    REQUIRED_P07_TOOLS,
    planned_command_templates,
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
        self.p07_bacterial = self._sequence_file("p07_bacterial.faa", ">p07-bac1 arbitrary_source_header\nMPEPTIDE\n>p07-review1\nMPEPTI\n")
        self.p07_archaeal = self._sequence_file("p07_archaeal.faa", ">p07-arc1\nMKKK\n")
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
                {"p07_sequence_id": "p07-arc1", "proteome_shard": "part_b", "target_id": "arc1", "source_proteome_path": "/machine/RS_GCF_000002.faa.gz", "target_length_from_p06": "4", "sequence_length": "4", "family_categories": "archaeal_patatin_like_pha_dep", "fasta_shard": str(self.p07_archaeal)},
                {"p07_sequence_id": "p07-bac1", "proteome_shard": "part_a", "target_id": "bac1", "source_proteome_path": "/machine/GB_GCF_000001.faa.gz", "target_length_from_p06": "8", "sequence_length": "8", "family_categories": "extracellular_pha_depolymerase_core", "fasta_shard": str(self.p07_bacterial)},
                {"p07_sequence_id": "p07-review1", "proteome_shard": "part_a", "target_id": "review1", "source_proteome_path": "/machine/GB_GCF_000001.faa.gz", "target_length_from_p06": "6", "sequence_length": "6", "family_categories": "extracellular_pha_depolymerase_core", "fasta_shard": str(self.p07_bacterial)},
            ],
        )
        write_tsv(
            self.p07_status,
            ("tool", "fasta_shard", "status"),
            [
                {"tool": tool, "fasta_shard": shard, "status": "completed"}
                for tool in REQUIRED_P07_TOOLS
                for shard in (str(self.p07_archaeal), str(self.p07_bacterial))
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
            ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path", "sequence_sha256", "evidence", "notes"),
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "model_sha256": "a" * 64, "seed_id": "seed-arc", "source_accession": "ARC1", "sequence_path": str(self.seed_archaeal), "sequence_sha256": self._sha256(self.seed_archaeal), "evidence": "E3", "notes": "archaeal seed"},
                {"family_category": "extracellular_pha_depolymerase_core", "model_sha256": "b" * 64, "seed_id": "seed-bac", "source_accession": "BAC1", "sequence_path": str(self.seed_bacterial), "sequence_sha256": self._sha256(self.seed_bacterial), "evidence": "E1", "notes": "bacterial seed"},
            ],
        )
        write_tsv(
            self.controls,
            ("family_category", "model_sha256", "control_id", "control_role", "sequence_path", "sequence_sha256", "evidence", "notes"),
            [
                {"family_category": "archaeal_patatin_like_pha_dep", "model_sha256": "a" * 64, "control_id": "control-arc", "control_role": "hard_negative", "sequence_path": str(self.control_archaeal), "sequence_sha256": self._sha256(self.control_archaeal), "evidence": "control", "notes": "archaeal control"},
                {"family_category": "extracellular_pha_depolymerase_core", "model_sha256": "b" * 64, "control_id": "control-bac", "control_role": "hard_negative", "sequence_path": str(self.control_bacterial), "sequence_sha256": self._sha256(self.control_bacterial), "evidence": "control", "notes": "bacterial control"},
            ],
        )
        self.taxonomy.write_text("GCF_000001\td__Bacteria;p__Test\nGCF_000002\td__Archaea;p__Test\n", encoding="utf-8")

    def _write_authoritative_core_reference_view(self) -> tuple[Path, Path]:
        """Replace synthetic core rows with the tracked legacy/core split schema."""
        source_families = (
            ("extracellular_mcl_pha_dep", "c" * 64, 6),
            ("extracellular_scl_pha_dep_type_I", "d" * 64, 6),
            ("extracellular_scl_pha_dep_type_II", "e" * 64, 5),
        )
        registry_rows = read_tsv(self.registry)
        registry_rows.extend(
            {"family_category": family, "approved_for_p06": "no", "scan_permission": "blocked", "model_sha256": model_sha}
            for family, model_sha, _ in source_families
        )
        write_tsv(self.registry, tuple(registry_rows[0].keys()), registry_rows)

        seed_rows = [row for row in read_tsv(self.seeds) if row["family_category"] != "extracellular_pha_depolymerase_core"]
        core_rows: list[dict[str, str]] = []
        index = 1
        for family, model_sha, count in source_families:
            for _ in range(count):
                accession = f"CORE{index:02d}"
                source = self._sequence_file(f"{accession}.faa", f">source-{accession}\nMPEPTIDE\n")
                seed_rows.append({
                    "family_category": family,
                    "model_sha256": model_sha,
                    "seed_id": f"legacy-{accession}",
                    "source_accession": accession,
                    "sequence_path": str(source),
                    "sequence_sha256": self._sha256(source),
                    "evidence": "E1",
                    "notes": "legacy extracellular source seed",
                })
                core_rows.append({
                    "family_category": "extracellular_pha_depolymerase_core",
                    "model_sha256": "b" * 64,
                    "seed_id": f"core-{accession}",
                    "source_accession": accession,
                    "sequence_path": str(source),
                })
                index += 1
        write_tsv(self.seeds, tuple(seed_rows[0].keys()), seed_rows)

        control_rows = [row for row in read_tsv(self.controls) if row["family_category"] != "extracellular_pha_depolymerase_core"]
        write_tsv(self.controls, tuple(control_rows[0].keys()), control_rows)
        core_seed_registry = self.root / "p05_extracellular_core_seed_registry.tsv"
        write_tsv(core_seed_registry, ("family_category", "model_sha256", "seed_id", "source_accession", "sequence_path"), core_rows)
        close_controls = self.root / "p05_extracellular_core_close_controls.tsv"
        residue_sha = hashlib.sha256(b"AAAA").hexdigest()
        write_tsv(
            close_controls,
            ("family_category", "source_accession", "sequence_path", "notes"),
            [{
                "family_category": "cutinase_like_non_target",
                "source_accession": "CLOSE1",
                "sequence_path": str(self.control_bacterial),
                "notes": f"Hard close non-target control; residue SHA256 is {residue_sha}.",
            }],
        )
        return core_seed_registry, close_controls

    def _add_route_family(self, family: str, candidate_count: int) -> None:
        fasta = self.root / f"{family}.faa"
        fasta.write_text("".join(f">p07-{family}-{index}\nM\n" for index in range(candidate_count)), encoding="utf-8")
        p06_rows = read_tsv(self.p06)
        p07_rows = read_tsv(self.p07_sequences)
        for index in range(candidate_count):
            target_id = f"{family}-{index}"
            p06_rows.append({
                "family_category": family,
                "proteome_shard": family,
                "target_id": target_id,
                "target_accession": target_id,
                "target_length": "1",
                "full_sequence_score": "200.0",
                "hmm_coverage": "0.9",
                "tier": "High-confidence",
            })
            p07_rows.append({
                "p07_sequence_id": f"p07-{family}-{index}",
                "proteome_shard": family,
                "target_id": target_id,
                "source_proteome_path": "/machine/GB_GCF_000001.faa.gz",
                "target_length_from_p06": "1",
                "sequence_length": "1",
                "family_categories": family,
                "fasta_shard": str(fasta),
            })
        write_tsv(self.p06, P06_FIELDS, p06_rows)
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, p07_rows)

        status_rows = read_tsv(self.p07_status)
        status_rows.extend({"tool": tool, "fasta_shard": str(fasta), "status": "completed"} for tool in REQUIRED_P07_TOOLS)
        write_tsv(self.p07_status, ("tool", "fasta_shard", "status"), status_rows)

        model_sha = hashlib.sha256(family.encode("ascii")).hexdigest()
        registry_rows = read_tsv(self.registry)
        registry_rows.append({"family_category": family, "approved_for_p06": "yes", "scan_permission": "approved", "model_sha256": model_sha})
        write_tsv(self.registry, tuple(registry_rows[0].keys()), registry_rows)
        seed = self._sequence_file(f"{family}.seed.faa", ">seed\nM\n")
        control = self._sequence_file(f"{family}.control.faa", ">control\nA\n")
        seed_rows = read_tsv(self.seeds)
        seed_rows.append({"family_category": family, "model_sha256": model_sha, "seed_id": f"seed-{family}", "source_accession": f"SEED-{family}", "sequence_path": str(seed), "sequence_sha256": self._sha256(seed), "evidence": "E1", "notes": "route seed"})
        write_tsv(self.seeds, tuple(seed_rows[0].keys()), seed_rows)
        control_rows = read_tsv(self.controls)
        control_rows.append({"family_category": family, "model_sha256": model_sha, "control_id": f"control-{family}", "control_role": "hard_negative", "sequence_path": str(control), "sequence_sha256": self._sha256(control), "evidence": "control", "notes": "route control"})
        write_tsv(self.controls, tuple(control_rows[0].keys()), control_rows)

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

    def test_p08_readme_exists_and_is_nonempty(self) -> None:
        readme = Path(__file__).resolve().parents[1] / "07_phylogeny" / "README.md"
        self.assertTrue(readme.is_file())
        self.assertTrue(readme.read_text(encoding="utf-8").strip())

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

    def test_cli_writes_planned_manifests_without_executing_phylogeny_tools(self) -> None:
        """The CLI combines Bac120/Ar53 taxonomy inputs but only plans P08 commands."""
        bac120_taxonomy = self.root / "bac120_taxonomy.tsv"
        ar53_taxonomy = self.root / "ar53_taxonomy.tsv"
        bac120_taxonomy.write_text("GCF_000001\td__Bacteria;p__Bac120\n", encoding="utf-8")
        ar53_taxonomy.write_text("GCF_000002\td__Archaea;p__Ar53\n", encoding="utf-8")
        bac120_tree = self.root / "bac120.tree"
        ar53_tree = self.root / "ar53.tree"
        bac120_tree.write_text("(GCF_000001:1);\n", encoding="utf-8")
        ar53_tree.write_text("(GCF_000002:1);\n", encoding="utf-8")
        cli_outdir = self.root / "cli-out"

        completed = subprocess.run(
            [
                sys.executable, "scripts/p08_prepare_phylogeny.py",
                "--candidate-table", str(self.p06),
                "--p07-sequence-manifest", str(self.p07_sequences),
                "--p07-status-table", str(self.p07_status),
                "--model-registry", str(self.registry),
                "--seed-registry", str(self.seeds),
                "--control-panel", str(self.controls),
                "--bac120-taxonomy", str(bac120_taxonomy),
                "--ar53-taxonomy", str(ar53_taxonomy),
                "--bac120-tree", str(bac120_tree),
                "--ar53-tree", str(ar53_tree),
                "--outdir", str(cli_outdir),
                "--include-tier", "High-confidence",
                "--include-tier", "Review",
                "--mafft-exe", "mafft-local",
                "--iqtree-exe", "iqtree2-local",
                "--fasttree-exe", "FastTree-local",
            ],
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        for name in (
            "p08_candidate_manifest.tsv", "p08_family_reference_manifest.tsv",
            "p08_preparation_summary.tsv", "p08_family_input_manifest.tsv",
            "p08_phylogeny_command_manifest.tsv", "p08_input_provenance.tsv",
        ):
            self.assertTrue((cli_outdir / "manifests" / name).is_file(), name)
        command_rows = read_tsv(cli_outdir / "manifests" / "p08_phylogeny_command_manifest.tsv")
        self.assertTrue(command_rows)
        self.assertTrue(all(row["command_status"] == "planned_not_run" for row in command_rows))
        self.assertTrue(all("mafft-local" in row["mafft_template"] for row in command_rows))
        self.assertFalse((cli_outdir / "alignments").exists())
        taxonomy_rows = read_tsv(cli_outdir / "gtdb_mapping" / "p08_taxonomy_join.tsv")
        self.assertEqual(
            {row["assembly_accession"]: row["taxonomy_lineage"] for row in taxonomy_rows},
            {"GCF_000001": "d__Bacteria;p__Bac120", "GCF_000002": "d__Archaea;p__Ar53"},
        )
        provenance_rows = read_tsv(cli_outdir / "manifests" / "p08_input_provenance.tsv")
        provenance = {row["input_role"]: row for row in provenance_rows}
        self.assertEqual(provenance["bac120_taxonomy"]["input_sha256"], self._sha256(bac120_taxonomy))
        self.assertEqual(provenance["ar53_taxonomy"]["input_sha256"], self._sha256(ar53_taxonomy))
        self.assertEqual(provenance["bac120_tree"]["input_sha256"], self._sha256(bac120_tree))
        self.assertEqual(provenance["ar53_tree"]["input_sha256"], self._sha256(ar53_tree))
        self.assertEqual(provenance["bac120_tree"]["input_usage"], "provenance_preflight_only_no_topology_read")
        self.assertEqual(provenance["ar53_tree"]["input_usage"], "provenance_preflight_only_no_topology_read")

        default_command = list(completed.args)
        include_index = default_command.index("--include-tier")
        del default_command[include_index:include_index + 4]
        default_outdir = self.root / "cli-default-out"
        default_command[default_command.index("--outdir") + 1] = str(default_outdir)
        default_completed = subprocess.run(
            default_command,
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(default_completed.returncode, 0, default_completed.stderr)
        default_candidates = read_tsv(default_outdir / "manifests" / "p08_candidate_manifest.tsv")
        self.assertNotIn("review1", [row["target_id"] for row in default_candidates])

    def test_cli_rejects_taxonomy_input_with_wrong_domain(self) -> None:
        bac120_taxonomy = self.root / "wrong_bac120_taxonomy.tsv"
        ar53_taxonomy = self.root / "ar53_taxonomy.tsv"
        bac120_taxonomy.write_text("GCF_000002\td__Archaea;p__WrongFile\n", encoding="utf-8")
        ar53_taxonomy.write_text("GCF_000001\td__Bacteria;p__WrongFile\n", encoding="utf-8")
        bac120_tree = self.root / "bac120.tree"
        ar53_tree = self.root / "ar53.tree"
        bac120_tree.write_text("(GCF_000001:1);\n", encoding="utf-8")
        ar53_tree.write_text("(GCF_000002:1);\n", encoding="utf-8")

        completed = subprocess.run(
            [
                sys.executable, "scripts/p08_prepare_phylogeny.py",
                "--candidate-table", str(self.p06),
                "--p07-sequence-manifest", str(self.p07_sequences),
                "--p07-status-table", str(self.p07_status),
                "--model-registry", str(self.registry),
                "--seed-registry", str(self.seeds),
                "--control-panel", str(self.controls),
                "--bac120-taxonomy", str(bac120_taxonomy),
                "--ar53-taxonomy", str(ar53_taxonomy),
                "--bac120-tree", str(bac120_tree),
                "--ar53-tree", str(ar53_tree),
                "--outdir", str(self.root / "wrong-domain-out"),
            ],
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Bac120 taxonomy", completed.stderr)

    def test_length_mismatch_blocks_and_writes_review(self) -> None:
        rows = read_tsv(self.p07_sequences)
        rows[0]["sequence_length"] = "5"
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, rows)
        with self.assertRaisesRegex(ValueError, "P06/P07 sequence length mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any("P06/P07 sequence length mismatch" in row["reason"] for row in blocks))

    def test_malformed_p06_or_p07_length_blocks_and_writes_review(self) -> None:
        p06_rows = read_tsv(self.p06)
        p07_rows = read_tsv(self.p07_sequences)
        p06_rows[0]["target_length"] = "not-a-length"
        p07_rows[0]["target_length_from_p06"] = "not-a-length"
        p07_rows[0]["sequence_length"] = "not-a-length"
        write_tsv(self.p06, P06_FIELDS, p06_rows)
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, p07_rows)
        with self.assertRaisesRegex(ValueError, "P06/P07 target length is not a positive integer"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "P06/P07 target length is not a positive integer" for row in blocks))

    def test_duplicate_p07_join_key_blocks_and_writes_review(self) -> None:
        rows = read_tsv(self.p07_sequences)
        duplicate = dict(rows[0])
        duplicate["p07_sequence_id"] = "p07-arc1-duplicate"
        rows.append(duplicate)
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, rows)
        with self.assertRaisesRegex(ValueError, "duplicate P07 join key"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "duplicate P07 join key" for row in blocks))

    def test_duplicate_p07_sequence_id_blocks_and_writes_review(self) -> None:
        rows = read_tsv(self.p07_sequences)
        duplicate = dict(rows[0])
        duplicate["proteome_shard"] = "another_shard"
        duplicate["target_id"] = "another_target"
        rows.append(duplicate)
        write_tsv(self.p07_sequences, P07_SEQUENCE_FIELDS, rows)
        with self.assertRaisesRegex(ValueError, "duplicate P07 sequence ID"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "duplicate P07 sequence ID" for row in blocks))

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
        rows = [row for row in read_tsv(self.p07_status) if not (row["tool"] == "InterProScan" and row["fasta_shard"] == str(self.p07_archaeal))]
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

    def test_seed_model_hash_mismatch_blocks_before_outputs_are_written(self) -> None:
        rows = read_tsv(self.seeds)
        rows[0]["model_sha256"] = "0" * 64
        write_tsv(self.seeds, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "P05 reference model SHA-256 mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        block = next(row for row in blocks if row["reason"] == "P05 reference model SHA-256 mismatch")
        self.assertEqual(block["family_category"], "archaeal_patatin_like_pha_dep")
        self.assertEqual(block["source_path"], str(self.seed_archaeal))
        self.assertIn("seed_id=seed-arc", block["notes"])
        self.assertFalse((self.outdir / "manifests" / "p08_family_input_manifest.tsv").exists())

    def test_control_model_hash_mismatch_blocks_before_outputs_are_written(self) -> None:
        rows = read_tsv(self.controls)
        rows[0]["model_sha256"] = "0" * 64
        write_tsv(self.controls, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "P05 reference model SHA-256 mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        block = next(row for row in blocks if row["reason"] == "P05 reference model SHA-256 mismatch")
        self.assertEqual(block["family_category"], "archaeal_patatin_like_pha_dep")
        self.assertEqual(block["source_path"], str(self.control_archaeal))
        self.assertIn("control_id=control-arc", block["notes"])
        self.assertFalse((self.outdir / "manifests" / "p08_phylogeny_command_manifest.tsv").exists())

    def test_authoritative_core_reference_view_derives_17_seeds_and_controls_without_subtype_output(self) -> None:
        core_seed_registry, _ = self._write_authoritative_core_reference_view()
        outputs = self._prepare()
        input_rows = read_tsv(outputs["family_input_manifest"])
        core_rows = [row for row in input_rows if row["family_category"] == "extracellular_pha_depolymerase_core"]
        core_seeds = [row for row in core_rows if row["record_kind"] == "seed"]
        core_controls = [row for row in core_rows if row["record_kind"] == "control"]
        self.assertEqual(len(core_seeds), 17)
        self.assertEqual(len(core_controls), 2)
        self.assertEqual({row["model_sha256"] for row in core_rows}, {"b" * 64})
        self.assertEqual({row["model_provenance"] for row in core_seeds}, {"derived_core_seed_registry"})
        self.assertEqual({row["model_provenance_source_path"] for row in core_seeds}, {str(core_seed_registry)})
        self.assertEqual({row["source_model_sha256"] for row in core_seeds}, {"c" * 64, "d" * 64, "e" * 64})
        command_families = {row["family_category"] for row in read_tsv(outputs["phylogeny_command_manifest"])}
        self.assertEqual(command_families, {"archaeal_patatin_like_pha_dep", "extracellular_pha_depolymerase_core"})

    def test_authoritative_core_seed_registry_hash_mismatch_blocks_and_writes_review(self) -> None:
        core_seed_registry, _ = self._write_authoritative_core_reference_view()
        rows = read_tsv(core_seed_registry)
        rows[0]["model_sha256"] = "0" * 64
        write_tsv(core_seed_registry, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "core seed registry model SHA-256 mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "core seed registry model SHA-256 mismatch" for row in blocks))

    def test_command_manifest_routes_199_200_and_2001_candidates_end_to_end(self) -> None:
        for family, count in (("route_199", 199), ("route_200", 200), ("route_2001", 2001)):
            self._add_route_family(family, count)
        rows = {row["family_category"]: row for row in read_tsv(self._prepare()["phylogeny_command_manifest"])}
        self.assertEqual(rows["route_199"]["candidate_input_record_count"], "199")
        self.assertEqual(rows["route_199"]["route"], "mafft_linsi_then_review")
        self.assertIn("mafft --localpair --maxiterate 1000", rows["route_199"]["mafft_template"])
        self.assertEqual(rows["route_200"]["candidate_input_record_count"], "200")
        self.assertEqual(rows["route_200"]["route"], "mafft_auto_then_review")
        self.assertIn("mafft --auto", rows["route_200"]["mafft_template"])
        self.assertEqual(rows["route_2001"]["candidate_input_record_count"], "2001")
        self.assertEqual(rows["route_2001"]["route"], "deterministic_representative_plan_then_fasttree_exploratory")
        self.assertIn("deterministic", rows["route_2001"]["representative_plan"])
        self.assertEqual(rows["route_2001"]["fasttree_template"], "FastTree -lg {representative_alignment_fasta} > {fasttree_tree}")
        self.assertTrue(all(rows[family]["command_status"] == "planned_not_run" for family in ("route_199", "route_200", "route_2001")))

    def test_route_family_size_boundaries(self) -> None:
        self.assertEqual(route_family_size(199), "mafft_linsi_then_review")
        self.assertEqual(route_family_size(200), "mafft_auto_then_review")
        self.assertEqual(route_family_size(2001), "deterministic_representative_plan_then_fasttree_exploratory")

    def test_planned_command_templates_route_at_required_boundaries(self) -> None:
        self.assertEqual(
            planned_command_templates(199)["mafft_template"],
            "mafft --localpair --maxiterate 1000 --thread {threads} --inputorder {input_fasta} > {alignment_fasta}",
        )
        self.assertEqual(
            planned_command_templates(200)["mafft_template"],
            "mafft --auto --thread {threads} --inputorder {input_fasta} > {alignment_fasta}",
        )
        large = planned_command_templates(2001)
        self.assertEqual(large["fasttree_template"], "FastTree -lg {representative_alignment_fasta} > {fasttree_tree}")
        self.assertIn("deterministic", large["representative_plan"])
        self.assertEqual(
            large["iqtree2_template"],
            "iqtree2 -s {alignment_fasta} -m TEST -B 1000 --prefix {iqtree_prefix}",
        )

    def test_preparation_writes_deterministic_family_fastas_and_planned_command_manifest(self) -> None:
        with patch("scripts.p08_prepare_phylogeny.subprocess.run") as run:
            outputs = self._prepare()
        run.assert_not_called()
        self.assertIn("family_input_manifest", outputs)
        self.assertIn("phylogeny_command_manifest", outputs)

        bacterial_fasta = self.outdir / "family_fastas" / "extracellular_pha_depolymerase_core.faa"
        self.assertEqual(
            bacterial_fasta.read_text(encoding="utf-8"),
            ">candidate|p07-bac1|extracellular_pha_depolymerase_core|GCF_000001\nMPEPTIDE\n"
            ">seed|seed-bac|extracellular_pha_depolymerase_core|BAC1\nMPEPTIDE\n"
            ">control|control-bac|extracellular_pha_depolymerase_core|control-bac\nAAAA\n",
        )
        input_rows = read_tsv(outputs["family_input_manifest"])
        bacterial_rows = [row for row in input_rows if row["family_category"] == "extracellular_pha_depolymerase_core"]
        self.assertEqual([row["record_kind"] for row in bacterial_rows], ["candidate", "seed", "control"])
        self.assertEqual([row["is_gtdb_candidate"] for row in bacterial_rows], ["yes", "no", "no"])
        self.assertTrue(all(row["input_fasta_path"] == str(bacterial_fasta) for row in bacterial_rows))
        self.assertEqual({row["input_sha256"] for row in bacterial_rows}, {self._sha256(bacterial_fasta)})
        self.assertTrue(all(row["evidence_boundary"] == "sequence_and_annotation_evidence_only_not_phenotype_proof" for row in input_rows))

        command_rows = read_tsv(outputs["phylogeny_command_manifest"])
        self.assertTrue(all(row["command_status"] == "planned_not_run" for row in command_rows))
        self.assertTrue(all(row["rooting_policy"] == "explicit_accessioned_outgroup_required; otherwise midpoint_display_only" for row in command_rows))
        bacterial_command = next(row for row in command_rows if row["family_category"] == "extracellular_pha_depolymerase_core")
        self.assertEqual(bacterial_command["candidate_input_record_count"], "1")
        self.assertEqual(bacterial_command["total_input_record_count"], "3")
        self.assertEqual(bacterial_command["route"], "mafft_linsi_then_review")
        self.assertEqual(
            bacterial_command["mafft_template"],
            "mafft --localpair --maxiterate 1000 --thread {threads} --inputorder {input_fasta} > {alignment_fasta}",
        )
        self.assertEqual(
            bacterial_command["iqtree2_template"],
            "iqtree2 -s {alignment_fasta} -m TEST -B 1000 --prefix {iqtree_prefix}",
        )
        self.assertEqual(bacterial_command["iqtree2_annotation"], "requires_independent_subset_and_outgroup_approval")

        summary_rows = read_tsv(outputs["preparation_summary"])
        bacterial_summary = next(row for row in summary_rows if row["family_category"] == "extracellular_pha_depolymerase_core")
        self.assertEqual(bacterial_summary["family_fasta_path"], str(bacterial_fasta))
        self.assertEqual(bacterial_summary["family_fasta_sha256"], self._sha256(bacterial_fasta))
        self.assertEqual(bacterial_summary["total_fasta_record_count"], "3")

    def test_missing_candidate_fasta_blocks_and_writes_review(self) -> None:
        self.p07_archaeal.unlink()
        with self.assertRaisesRegex(ValueError, "candidate FASTA unreadable"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "candidate FASTA unreadable" for row in blocks))

    def test_absent_p07_sequence_id_blocks_and_writes_review(self) -> None:
        self.p07_archaeal.write_text(">another-record\nMKKK\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "candidate FASTA missing p07_sequence_id"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "candidate FASTA missing p07_sequence_id" for row in blocks))

    def test_candidate_fasta_length_mismatch_blocks_and_writes_review(self) -> None:
        self.p07_archaeal.write_text(">p07-arc1\nMKKKK\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "candidate FASTA sequence length mismatch"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "candidate FASTA sequence length mismatch" for row in blocks))

    def test_malformed_or_checksum_conflicting_reference_blocks_and_writes_review(self) -> None:
        self.seed_archaeal.write_text("MKKK\n", encoding="utf-8")
        rows = read_tsv(self.seeds)
        rows[0]["sequence_sha256"] = self._sha256(self.seed_archaeal)
        write_tsv(self.seeds, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "reference FASTA malformed"):
            self._prepare()
        blocks = read_tsv(self.outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "reference FASTA malformed" for row in blocks))

        checksum_outdir = self.root / "checksum-conflict"
        self.seed_archaeal.write_text(">seed-a\nMKKK\n", encoding="utf-8")
        self._write_complete_fixture()
        rows = read_tsv(self.seeds)
        rows[0]["sequence_sha256"] = "0" * 64
        write_tsv(self.seeds, tuple(rows[0].keys()), rows)
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self._prepare(outdir=checksum_outdir)
        checksum_blocks = read_tsv(checksum_outdir / "review" / "p08_blocked_records.tsv")
        self.assertTrue(any(row["reason"] == "SHA-256 mismatch" for row in checksum_blocks))

    def test_prepare_p08_inputs_never_calls_subprocess_run(self) -> None:
        with patch("scripts.p08_prepare_phylogeny.subprocess.run") as run:
            self._prepare()
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
