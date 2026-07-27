from __future__ import annotations

import csv
import tempfile
import textwrap
from pathlib import Path
import unittest

from scripts import p05_plan_family_profiles as p05
from scripts.p05_plan_family_profiles import build_family_profile_plan, plan_family_profiles, summarize_family_profile_plan


class P05PlanFamilyProfilesTest(unittest.TestCase):
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

    def test_plan_family_profiles_separates_hmm_ready_and_anchor_sets(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(seed_id="seed-a3", family_category="phaZ7_like", source_accession="A3", evidence_level="E2"),
            self._seed_row(seed_id="seed-b1", family_category="phaZd_like", source_accession="B1"),
            self._seed_row(seed_id="seed-b2", family_category="phaZd_like", source_accession="B2", evidence_level="E2"),
        ]

        plan_rows = plan_family_profiles(rows, minimum_independent_seeds=3)
        summary_rows = summarize_family_profile_plan(plan_rows)

        self.assertEqual(len(plan_rows), 2)
        self.assertEqual({row["profile_strategy"] for row in plan_rows}, {"build_hmm", "anchor_set"})
        hmm_ready = next(row for row in plan_rows if row["family_category"] == "phaZ7_like")
        anchor = next(row for row in plan_rows if row["family_category"] == "phaZd_like")
        self.assertEqual(hmm_ready["eligible_for_hmm"], "yes")
        self.assertEqual(hmm_ready["independent_qualifying_accession_count"], "3")
        self.assertEqual(anchor["eligible_for_hmm"], "no")
        self.assertEqual(anchor["independent_qualifying_accession_count"], "2")
        self.assertIn({"kind": "profile_strategy", "name": "build_hmm", "count": "1"}, summary_rows)
        self.assertIn({"kind": "profile_strategy", "name": "anchor_set", "count": "1"}, summary_rows)

    def test_plan_family_profiles_counts_archaeal_e3_annotation_seeds(self) -> None:
        rows = [
            self._seed_row(
                seed_id="arc-a1",
                reference_library="archaea_literature_supported",
                taxonomic_domain="Archaea",
                family_category="archaeal_patatin_like_pha_dep",
                seed_name="Annot A1",
                evidence_level="E3",
                source_database="NCBI Protein",
                source_accession="AHZ23723.1",
                organism="Haloferax mediterranei ATCC 33500",
                taxon_id="523841",
                sequence_length_aa="442",
                sequence_path="seeds/archaea/AHZ23723.1.faa",
            ),
            self._seed_row(
                seed_id="arc-a2",
                reference_library="archaea_literature_supported",
                taxonomic_domain="Archaea",
                family_category="archaeal_patatin_like_pha_dep",
                seed_name="Annot A2",
                evidence_level="E3",
                source_database="NCBI Protein",
                source_accession="AHB64615.1",
                organism="Haloarcula hispanica N601",
                taxon_id="1417673",
                sequence_length_aa="474",
                sequence_path="seeds/archaea/AHB64615.1.faa",
            ),
            self._seed_row(
                seed_id="arc-a3",
                reference_library="archaea_literature_supported",
                taxonomic_domain="Archaea",
                family_category="archaeal_patatin_like_pha_dep",
                seed_name="Annot A3",
                evidence_level="E3",
                source_database="NCBI Protein",
                source_accession="CCQ36014.1",
                organism="Natronomonas moolapensis 8.8.11",
                taxon_id="268739",
                sequence_length_aa="323",
                sequence_path="seeds/archaea/CCQ36014.1.faa",
            ),
        ]

        plan_rows = plan_family_profiles(rows, minimum_independent_seeds=3)

        self.assertEqual(len(plan_rows), 1)
        self.assertEqual(plan_rows[0]["profile_strategy"], "build_hmm")
        self.assertEqual(plan_rows[0]["independent_qualifying_accession_count"], "3")
        self.assertEqual(plan_rows[0]["qualifying_source_accessions"], "AHB64615.1;AHZ23723.1;CCQ36014.1")

    def test_plan_excludes_boundary_candidates_from_hmm_seed_count(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a1", family_category="phaZ7_like", source_accession="A1"),
            self._seed_row(seed_id="seed-a2", family_category="phaZ7_like", source_accession="A2"),
            self._seed_row(
                seed_id="seed-boundary",
                family_category="phaZ7_like",
                source_accession="A3",
                evidence_level="E2",
                profile_seed_status="boundary_candidate",
            ),
        ]

        plan_rows = plan_family_profiles(rows, minimum_independent_seeds=3)

        self.assertEqual(plan_rows[0]["profile_strategy"], "anchor_set")
        self.assertEqual(plan_rows[0]["qualifying_seed_row_count"], "2")
        self.assertEqual(plan_rows[0]["qualifying_source_accessions"], "A1;A2")

    def test_build_family_profile_plan_filters_to_keep_now_families(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "reference_library.normalized.tsv"
            classification = root / "p05_family_keep_now.tsv"
            manifest.write_text(
                textwrap.dedent(
                    """
                    seed_id\treference_library\ttaxonomic_domain\tfamily_category\tseed_name\tevidence_level\tsource_database\tsource_accession\torganism\ttaxon_id\tretrieval_date\tsequence_format\tsequence_length_aa\tsequence_path
                    seed-a1\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed A1\tE1\tUniProtKB\tA1\tOrg A1\t1\t2026-07-24\tfaa\t101\tseeds/a1.faa
                    seed-a2\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed A2\tE1\tUniProtKB\tA2\tOrg A2\t2\t2026-07-24\tfaa\t102\tseeds/a2.faa
                    seed-a3\tbacteria_high_confidence\tBacteria\tphaZ7_like\tSeed A3\tE2\tUniProtKB\tA3\tOrg A3\t3\t2026-07-24\tfaa\t103\tseeds/a3.faa
                    seed-b1\tbacteria_high_confidence\tBacteria\tphaZd_like\tSeed B1\tE1\tUniProtKB\tB1\tOrg B1\t4\t2026-07-24\tfaa\t104\tseeds/b1.faa
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )
            p05.write_tsv(
                classification,
                [
                    {
                        "seed_id": "seed-a1",
                        "family_category": "phaZ7_like",
                        "source_accession": "A1",
                        "priority_status": "keep_now",
                        "reason": "active main class",
                    }
                ],
                ("seed_id", "family_category", "source_accession", "priority_status", "reason"),
            )

            outputs = build_family_profile_plan(
                manifest,
                root / "manifests",
                classification_path=classification,
            )

            for path in outputs.values():
                self.assertTrue(path.is_file())

            with outputs["plan"].open("r", encoding="utf-8", newline="") as handle:
                plan_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["hmm_ready"].open("r", encoding="utf-8", newline="") as handle:
                hmm_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["anchor_sets"].open("r", encoding="utf-8", newline="") as handle:
                anchor_rows = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["summary"].open("r", encoding="utf-8", newline="") as handle:
                summary_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(plan_rows), 1)
            self.assertEqual(len(hmm_rows), 1)
            self.assertEqual(len(anchor_rows), 0)
            self.assertEqual(hmm_rows[0]["family_category"], "phaZ7_like")
            self.assertIn({"kind": "total", "name": "families", "count": "1"}, summary_rows)

    def test_build_family_profile_plan_rejects_missing_keep_now_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "reference_library.normalized.tsv"
            classification = root / "p05_family_keep_now.tsv"
            p05.write_tsv(
                manifest,
                [self._seed_row(family_category="phaZ7_like", source_accession="A1")],
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
                        "family_category": "not_in_manifest",
                        "priority_status": "keep_now",
                    }
                ],
                ("family_category", "priority_status"),
            )

            with self.assertRaisesRegex(ValueError, "not present in the reference manifest"):
                build_family_profile_plan(
                    manifest,
                    root / "manifests",
                    classification_path=classification,
                )

    def test_plan_family_profiles_rejects_mixed_domains(self) -> None:
        rows = [
            self._seed_row(seed_id="seed-a", family_category="phaZ7_like", taxonomic_domain="Bacteria", source_accession="A"),
            self._seed_row(seed_id="seed-b", family_category="phaZ7_like", taxonomic_domain="Archaea", source_accession="B"),
        ]

        with self.assertRaisesRegex(ValueError, "spans multiple taxonomic domains"):
            plan_family_profiles(rows)


if __name__ == "__main__":
    unittest.main()
