from __future__ import annotations

import csv
import tempfile
from pathlib import Path
import unittest

from scripts import p05_catalog_hmm_models as catalog


class P05CatalogHmmModelsTest(unittest.TestCase):
    def _write_tsv(self, path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _make_fixture(self, root: Path) -> tuple[Path, Path, Path, Path, Path, Path]:
        reference_fasta = root / "A1.faa"
        reference_fasta.write_text(">seed-a1|A1\nMST\n", encoding="ascii")
        reference = root / "reference.tsv"
        self._write_tsv(
            reference,
            ["family_category", "seed_id", "source_accession", "sequence_path", "organism"],
            [{"family_category": "test_family", "seed_id": "seed-a1", "source_accession": "A1", "sequence_path": reference_fasta.as_posix(), "organism": "Test organism"}],
        )
        commands = root / "commands.tsv"
        self._write_tsv(commands, ["family_category", "hmmbuild_command"], [{"family_category": "test_family", "hmmbuild_command": "hmmbuild --amino test.hmm test.aligned.faa"}])
        models = root / "hmms"
        bundles = root / "bundles"
        alignments = root / "alignments"
        models.mkdir()
        bundles.mkdir()
        alignments.mkdir()
        (models / "test_family.hmm").write_text(
            "HMMER3/f [3.4 | Aug 2023]\nNAME  test_family.aligned\nLENG  3\nDATE  Sun Jul 26 19:37:10 2026\nNSEQ  1\nEFFN  1.000000\n//\nHMM\n",
            encoding="ascii",
        )
        (bundles / "test_family.faa").write_text(">seed-a1|A1\nMST\n", encoding="ascii")
        (alignments / "test_family.aligned.faa").write_text(">seed-a1|A1\nMST\n", encoding="ascii")
        return reference, commands, models, bundles, alignments, root / "out"

    def test_catalog_writes_checksum_locked_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            reference, commands, models, bundles, alignments, outdir = self._make_fixture(Path(tmp))
            outputs = catalog.catalog_hmm_models(reference, commands, models, bundles, alignments, outdir)

            with outputs["model_registry"].open("r", encoding="utf-8", newline="") as handle:
                model_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(model_row["family_category"], "test_family")
            self.assertEqual(model_row["hmmer_version"], "3.4")
            self.assertEqual(model_row["seed_sequence_count"], "1")
            self.assertEqual(model_row["approved_for_p06"], "no")
            self.assertEqual(len(model_row["model_sha256"]), 64)

            with outputs["seed_registry"].open("r", encoding="utf-8", newline="") as handle:
                seed_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(seed_row["source_accession"], "A1")
            self.assertEqual(seed_row["bundle_membership_verified"], "yes")
            self.assertEqual(len(seed_row["sequence_sha256"]), 64)

    def test_catalog_rejects_nseq_bundle_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            reference, commands, models, bundles, alignments, outdir = self._make_fixture(Path(tmp))
            (models / "test_family.hmm").write_text(
                "HMMER3/f [3.4 | Aug 2023]\nNAME  test_family.aligned\nLENG  3\nDATE  Sun Jul 26 19:37:10 2026\nNSEQ  2\nEFFN  1.000000\n//\nHMM\n",
                encoding="ascii",
            )
            with self.assertRaisesRegex(ValueError, "NSEQ"):
                catalog.catalog_hmm_models(reference, commands, models, bundles, alignments, outdir)


if __name__ == "__main__":
    unittest.main()
