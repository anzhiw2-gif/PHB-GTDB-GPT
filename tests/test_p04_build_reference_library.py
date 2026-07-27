from __future__ import annotations

import csv
import tempfile
import textwrap
from pathlib import Path
import unittest

from scripts.p04_build_reference_library import (
    CANONICAL_EVIDENCE_LEVELS,
    CANONICAL_FAMILY_CATEGORIES,
    build_reference_library,
    load_reference_manifest,
    normalize_rows,
    summarize_rows,
)


class P04BuildReferenceLibraryTest(unittest.TestCase):
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

    def test_load_reference_manifest_reads_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                textwrap.dedent(
                    """
                    seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path
                    seed-1\tbacteria_high_confidence\tBacteria\tintracellular_phaZ_no_lipase_box\tPhaZ alpha\tE1\tUniProtKB\tP00001.1\tBacillus sp.\t12345\t2026-07-24\tfaa\t312\tseeds/phaZ_alpha.faa
                    seed-2\tbacteria_high_confidence\tBacteria\textracellular_scl_pha_dep_type_I\tScl-PHA depolymerase\tE2\tNCBI Protein\tNP_000002.1\tPseudomonas sp.\t54321\t2026-07-24\tfasta\t280\tseeds/scl_type_I.faa
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            rows = load_reference_manifest(manifest)

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["family_category"], "intracellular_phaZ_no_lipase_box")
            self.assertEqual(rows[1]["evidence_level"], "E2")
            self.assertEqual(rows[0]["taxon_id"], "12345")
            self.assertEqual(rows[1]["sequence_format"], "fasta")

    def test_normalize_rows_orders_by_family_evidence_and_seed(self) -> None:
        rows = [
            {
                "seed_id": "seed-b",
                "reference_library": "bacteria_high_confidence",
                "taxonomic_domain": "Bacteria",
                "family_category": "extracellular_scl_pha_dep_type_II",
                "seed_name": "b",
                "evidence_level": "E2",
                "source_database": "UniProtKB",
                "source_accession": "B",
                "organism": "Org B",
                "taxon_id": "2",
                "retrieval_date": "2026-07-24",
                "sequence_format": "faa",
                "sequence_length_aa": "120",
                "sequence_path": "b.faa",
            },
            {
                "seed_id": "seed-a",
                "reference_library": "bacteria_high_confidence",
                "taxonomic_domain": "Bacteria",
                "family_category": "extracellular_scl_pha_dep_type_I",
                "seed_name": "a",
                "evidence_level": "E1",
                "source_database": "NCBI Protein",
                "source_accession": "A",
                "organism": "Org A",
                "taxon_id": "1",
                "retrieval_date": "2026-07-24",
                "sequence_format": "faa",
                "sequence_length_aa": "121",
                "sequence_path": "a.faa",
            },
        ]

        normalized = normalize_rows(rows)

        self.assertEqual([row["seed_id"] for row in normalized], ["seed-a", "seed-b"])

    def test_summarize_rows_counts_each_family_and_evidence_level(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-1", family_category=CANONICAL_FAMILY_CATEGORIES[0], evidence_level=CANONICAL_EVIDENCE_LEVELS[0], source_accession="A"),
            self._seed_row(seed_id="seed-2", family_category=CANONICAL_FAMILY_CATEGORIES[0], evidence_level=CANONICAL_EVIDENCE_LEVELS[1], source_accession="B"),
        ]

        summary = summarize_rows(rows)

        self.assertIn({"kind": "family", "name": CANONICAL_FAMILY_CATEGORIES[0], "count": 2}, summary)
        self.assertIn({"kind": "reference_library", "name": "bacteria_high_confidence", "count": 2}, summary)
        self.assertIn({"kind": "taxonomic_domain", "name": "Bacteria", "count": 2}, summary)
        self.assertIn({"kind": "evidence", "name": "E1", "count": 1}, summary)
        self.assertIn({"kind": "evidence", "name": "E2", "count": 1}, summary)

    def test_build_reference_library_writes_normalized_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.tsv"
            manifest.write_text(
                textwrap.dedent(
                    """
                    seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\tnotes
                    seed-2\tbacteria_high_confidence\tBacteria\tphaZd_like\tSeed two\tE2\tUniProtKB\tQ9XYZ2.1\tOrg 2\t2002\t2026-07-24\tfaa\t245\tseeds/two.faa\tsecond
                    seed-1\tbacteria_high_confidence\tBacteria\tphaZd_like\tSeed one\tE1\tNCBI Protein\tABC111.1\tOrg 1\t2001\t2026-07-24\tfaa\t244\tseeds/one.faa\tfirst
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            outputs = build_reference_library(manifest, root / "out")

            self.assertTrue(outputs["normalized"].is_file())
            self.assertTrue(outputs["bacteria"].is_file())
            self.assertTrue(outputs["archaea"].is_file())
            self.assertTrue(outputs["summary"].is_file())
            with outputs["normalized"].open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual([row["seed_id"] for row in rows], ["seed-1", "seed-2"])
            self.assertEqual(rows[0]["source_accession"], "ABC111.1")
            self.assertEqual(rows[1]["notes"], "second")

    def test_load_reference_manifest_rejects_bad_category(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tunknown_family\tSeed\tE1\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "invalid family_category"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_bad_evidence_level(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed\tE9\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "invalid evidence_level"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_bad_source_database(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed\tE1\tPubMed\tPMID:123\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "invalid source_database"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_excluded_without_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed\tExcluded\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "exclusion_reason is required"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_accepts_literature_supported_archaeal_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\tpmid\tdoi\tliterature_support_scope\n"
                "archaea-1\tarchaea_literature_supported\tArchaea\tarchaeal_patatin_like_pha_dep\tPhaZh1\tE1\tNCBI Protein\tAFK21580.1\tHaloferax mediterranei ATCC 33500\t523841\t2026-07-25\tfaa\t321\tseeds/archaea/AFK21580.1.faa\t25710370\t10.1128/AEM.04269-14\tdirect native PHA depolymerase assay\n",
                encoding="utf-8",
            )

            rows = load_reference_manifest(manifest)

            self.assertEqual(rows[0]["reference_library"], "archaea_literature_supported")
            self.assertEqual(rows[0]["taxonomic_domain"], "Archaea")

    def test_load_reference_manifest_accepts_archaeal_annotation_seed_without_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\tliterature_support_scope\n"
                "archaea-2\tarchaea_literature_supported\tArchaea\tarchaeal_patatin_like_pha_dep\tAnnot. depolymerase\tE3\tNCBI Protein\tCCQ36014.1\tNatronomonas moolapensis 8.8.11\t268739\t2026-07-26\tfaa\t323\tseeds/archaea/CCQ36014.1.faa\tNCBI Protein annotation only; no direct literature support\n",
                encoding="utf-8",
            )

            rows = load_reference_manifest(manifest)

            self.assertEqual(rows[0]["evidence_level"], "E3")

    def test_load_reference_manifest_rejects_unknown_profile_seed_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\tprofile_seed_status\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed\tE1\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\tnot_a_valid_status\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "invalid profile_seed_status"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_archaeal_seed_without_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "archaea-1\tarchaea_literature_supported\tArchaea\tarchaeal_patatin_like_pha_dep\tPhaZh1\tE2\tNCBI Protein\tAFK21580.1\tHaloferax mediterranei ATCC 33500\t523841\t2026-07-25\tfaa\t321\tseeds/archaea/AFK21580.1.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "require PMID, DOI, or PMCID"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_bacterial_low_confidence_seed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed\tE3\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "bacteria_high_confidence rows require E1 or E2"):
                load_reference_manifest(manifest)

    def test_load_reference_manifest_rejects_domain_library_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path\n"
                "seed-1\tarchaea_literature_supported\tBacteria\tphaZ7_like\tSeed\tE1\tUniProtKB\tABC123.1\tOrg\t123\t2026-07-24\tfaa\t100\tseed.faa\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "requires reference_library"):
                load_reference_manifest(manifest)


if __name__ == "__main__":
    unittest.main()
