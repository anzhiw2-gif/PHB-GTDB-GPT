from __future__ import annotations

import csv
import tempfile
import textwrap
from pathlib import Path
import unittest

from scripts import p05_plan_family_profiles as p05


class P05BuildFamilyProfilesTest(unittest.TestCase):
    def _seed_row(self, **overrides: str) -> dict[str, str]:
        row = {
            "seed_id": "seed-1",
            "reference_library": "bacteria_high_confidence",
            "taxonomic_domain": "Bacteria",
            "family_category": "intracellular_phaZ_no_lipase_box",
            "seed_name": "PhaZ alpha",
            "evidence_level": "E1",
            "source_database": "UniProtKB",
            "source_accession": "P00001.1",
            "organism": "Bacillus sp.",
            "taxon_id": "12345",
            "retrieval_date": "2026-07-24",
            "sequence_format": "faa",
            "sequence_length_aa": "312",
            "sequence_path": "seeds/phaZ_alpha.faa",
            "notes": "",
        }
        row.update(overrides)
        return row

    def _write_manifest_and_plan(self, root: Path, rows: list[dict[str, str]]) -> tuple[Path, Path]:
        manifest = root / "reference_library.normalized.tsv"
        plan = root / "p05_family_profile_plan.tsv"
        for index, row in enumerate(rows, start=1):
            sequence_path = Path(row["sequence_path"])
            sequence_file = sequence_path if sequence_path.is_absolute() else root / sequence_path
            sequence_file.parent.mkdir(parents=True, exist_ok=True)
            sequence_file.write_text(
                f">input-{index}\n{'M' * (10 + index)}\n",
                encoding="utf-8",
            )
        p05.write_tsv(
            manifest,
            rows,
            (
                "seed_id",
                "reference_library",
                "taxonomic_domain",
                "family_category",
                "seed_name",
                "evidence_level",
                "source_database",
                "source_accession",
                "organism",
                "taxon_id",
                "retrieval_date",
                "sequence_format",
                "sequence_length_aa",
                "sequence_path",
                "notes",
            ),
        )
        plan_rows = p05.plan_family_profiles(rows, minimum_independent_seeds=3)
        p05.write_tsv(plan, plan_rows, p05.PLAN_FIELDNAMES)
        return manifest, plan

    def test_build_family_profile_scaffold_writes_empty_queue_when_no_family_qualifies(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZd_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZd_like", source_accession="A2"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)

            outputs = p05.build_family_profile_scaffold(manifest, plan, root / "build_manifests")

            for path in outputs.values():
                self.assertTrue(path.is_file())

            with outputs["queue"].open("r", encoding="utf-8", newline="") as handle:
                queue_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(queue_rows, [])
            self.assertIn({"kind": "total", "name": "families_in_plan", "count": "1"}, summary_rows)
            self.assertIn({"kind": "total", "name": "eligible_families", "count": "0"}, summary_rows)
            self.assertIn({"kind": "total", "name": "build_queue_rows", "count": "0"}, summary_rows)

    def test_main_build_scaffold_mode_writes_outputs(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3", evidence_level="E2"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)
            outdir = root / "build_manifests"

            exit_code = p05.main(
                [
                    "--manifest",
                    str(manifest),
                    "--plan-path",
                    str(plan),
                    "--outdir",
                    str(outdir),
                    "--build-scaffold",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((outdir / "p05_family_hmm_build_scaffold_queue.tsv").is_file())
            self.assertTrue((outdir / "p05_family_hmm_build_scaffold_summary.tsv").is_file())

    def test_main_build_scaffold_mode_generates_plan_when_missing(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3", evidence_level="E2"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "reference_library.normalized.tsv"
            classification = root / "p05_family_keep_now.tsv"
            for index, row in enumerate(rows, start=1):
                sequence_path = Path(row["sequence_path"])
                sequence_file = sequence_path if sequence_path.is_absolute() else root / sequence_path
                sequence_file.parent.mkdir(parents=True, exist_ok=True)
                sequence_file.write_text(
                    f">input-{index}\n{'M' * (10 + index)}\n",
                    encoding="utf-8",
                )
            p05.write_tsv(
                manifest,
                rows,
                (
                    "seed_id",
                    "reference_library",
                    "taxonomic_domain",
                    "family_category",
                    "seed_name",
                    "evidence_level",
                    "source_database",
                    "source_accession",
                    "organism",
                    "taxon_id",
                    "retrieval_date",
                    "sequence_format",
                    "sequence_length_aa",
                    "sequence_path",
                    "notes",
                ),
            )
            p05.write_tsv(
                classification,
                [
                    {
                        "family_category": "phaZ7_like",
                        "priority_status": "keep_now",
                    }
                ],
                ("family_category", "priority_status"),
            )

            outdir = root / "build_manifests"
            exit_code = p05.main(
                [
                    "--manifest",
                    str(manifest),
                    "--outdir",
                    str(outdir),
                    "--plan-path",
                    str(root / "missing_plan.tsv"),
                    "--family-classification",
                    str(classification),
                    "--build-scaffold",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((outdir / "p05_family_profile_plan.tsv").is_file())
            self.assertTrue((outdir / "p05_family_hmm_build_scaffold_queue.tsv").is_file())
            self.assertTrue((outdir / "p05_family_hmm_build_scaffold_summary.tsv").is_file())

    def test_build_family_profile_scaffold_queues_mafft_hmmer_family(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3", evidence_level="E2"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)

            outputs = p05.build_family_profile_scaffold(manifest, plan, root / "build_manifests")

            with outputs["queue"].open("r", encoding="utf-8", newline="") as handle:
                queue_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(queue_rows), 1)
            self.assertEqual(queue_rows[0]["family_category"], "phaZ7_like")
            self.assertEqual(queue_rows[0]["independent_qualifying_accession_count"], "3")
            self.assertEqual(queue_rows[0]["alignment_tool"], "MAFFT")
            self.assertEqual(queue_rows[0]["alignment_mode"], "L-INS-i")
            self.assertEqual(queue_rows[0]["hmm_build_tool"], "hmmbuild")
            self.assertEqual(queue_rows[0]["calibration_search_tool"], "hmmsearch")
            self.assertEqual(queue_rows[0]["bundled_sequence_count"], "3")
            bundle_path = Path(queue_rows[0]["seed_bundle_path"])
            self.assertTrue(bundle_path.is_file())
            self.assertTrue(bundle_path.read_text(encoding="utf-8").startswith(">seed-a1|A1\n"))
            self.assertIn({"kind": "total", "name": "eligible_families", "count": "1"}, summary_rows)
            self.assertIn({"kind": "total", "name": "build_queue_rows", "count": "1"}, summary_rows)

    def test_build_family_profile_scaffold_materializes_deterministic_bundle(self) -> None:
        rows = [
            self._seed_row(
                seed_id="seed-z2",
                family_category="phaZ7_like",
                source_accession="Z2",
                sequence_path="seeds/z2.faa",
            ),
            self._seed_row(
                seed_id="seed-a1",
                family_category="phaZ7_like",
                source_accession="A1",
                sequence_path="seeds/a1.faa",
            ),
            self._seed_row(
                seed_id="seed-m3",
                family_category="phaZ7_like",
                source_accession="M3",
                sequence_path="seeds/m3.faa",
            ),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)
            bundle_dir = root / "seed_bundles"

            outputs = p05.build_family_profile_scaffold(
                manifest,
                plan,
                root / "build_manifests",
                bundle_dir=bundle_dir,
            )

            with outputs["queue"].open("r", encoding="utf-8", newline="") as handle:
                queue_row = next(csv.DictReader(handle, delimiter="\t"))
            bundle_path = bundle_dir / "phaZ7_like.faa"
            self.assertEqual(queue_row["seed_bundle_path"], bundle_path.as_posix())
            self.assertEqual(queue_row["bundled_sequence_count"], "3")
            self.assertEqual(
                bundle_path.read_text(encoding="utf-8"),
                f">seed-a1|A1\n{'M' * 12}\n>seed-m3|M3\n{'M' * 13}\n>seed-z2|Z2\n{'M' * 11}\n",
            )

    def test_build_family_profile_scaffold_resolves_repo_root_relative_paths(self) -> None:
        rows = [
            self._seed_row(
                seed_id="seed-a1",
                family_category="phaZ7_like",
                source_accession="A1",
                sequence_path="01_reference_library/seeds/a1.faa",
            ),
            self._seed_row(
                seed_id="seed-a2",
                family_category="phaZ7_like",
                source_accession="A2",
                sequence_path="01_reference_library/seeds/a2.faa",
            ),
            self._seed_row(
                seed_id="seed-a3",
                family_category="phaZ7_like",
                source_accession="A3",
                sequence_path="01_reference_library/seeds/a3.faa",
            ),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)
            bundle_dir = root / "seed_bundles"

            outputs = p05.build_family_profile_scaffold(
                manifest,
                plan,
                root / "build_manifests",
                bundle_dir=bundle_dir,
            )

            bundle_path = bundle_dir / "phaZ7_like.faa"
            self.assertTrue(bundle_path.is_file())
            self.assertEqual(
                bundle_path.read_text(encoding="utf-8"),
                f">seed-a1|A1\n{'M' * 11}\n>seed-a2|A2\n{'M' * 12}\n>seed-a3|A3\n{'M' * 13}\n",
            )
            self.assertTrue(outputs["queue"].is_file())

    def test_build_family_profile_scaffold_keeps_bundle_path_consistent_on_rerun(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)
            outdir = root / "build_manifests"
            bundle_dir = root / "seed_bundles"

            first = p05.build_family_profile_scaffold(manifest, plan, outdir, bundle_dir=bundle_dir)
            with first["queue"].open("r", encoding="utf-8", newline="") as handle:
                first_queue = list(csv.DictReader(handle, delimiter="\t"))
            first_content = Path(first_queue[0]["seed_bundle_path"]).read_text(encoding="utf-8")
            second = p05.build_family_profile_scaffold(manifest, plan, outdir, bundle_dir=bundle_dir)
            with second["queue"].open("r", encoding="utf-8", newline="") as handle:
                second_queue = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(first_queue[0]["seed_bundle_path"], (bundle_dir / "phaZ7_like.faa").as_posix())
            self.assertEqual(second_queue[0]["seed_bundle_path"], first_queue[0]["seed_bundle_path"])
            self.assertEqual(Path(second_queue[0]["seed_bundle_path"]).read_text(encoding="utf-8"), first_content)

    def test_build_family_profile_scaffold_rejects_plan_manifest_mismatch(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3"),
        ]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, plan = self._write_manifest_and_plan(root, rows)

            manifest.write_text(
                textwrap.dedent(
                    """
                    seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\tnotes
                    seed-a1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed A1\tE1\tUniProtKB\tA1\tOrg A1\t1\t2026-07-24\tfaa\t101\tseeds/a1.faa\t
                    seed-a2\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed A2\tE1\tUniProtKB\tA2\tOrg A2\t2\t2026-07-24\tfaa\t102\tseeds/a2.faa\t
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "plan and manifest disagree"):
                p05.build_family_profile_scaffold(manifest, plan, root / "build_manifests")


if __name__ == "__main__":
    unittest.main()
