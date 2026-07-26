from __future__ import annotations

import gzip
import tempfile
import textwrap
from pathlib import Path
import unittest

from scripts.p03_predict_proteomes import (
    _gene_translation,
    discover_genome_files,
    infer_accession_from_genome_path,
    load_prediction_policy,
    output_paths_for_genome,
    predict_genome,
    run_p03_prediction,
    write_prediction_outputs,
)


class P03PredictProteomesTest(unittest.TestCase):
    def test_gene_translation_uses_pyrodigal_translate_method(self) -> None:
        class PyrodigalLikeGene:
            def translate(self, include_stop: bool = True) -> str:
                return "MKK*" if include_stop else "MKK"

        self.assertEqual(_gene_translation(PyrodigalLikeGene()), "MKK")

    def test_infer_accession_uses_real_gtdb_filename(self) -> None:
        path = Path("00_raw_gtdb_r232/genomes/GCF/033/239/625/GCF_033239625.1_genomic.fna.gz")

        self.assertEqual(infer_accession_from_genome_path(path), "GCF_033239625.1")

    def test_discover_genome_files_finds_sorted_fasta_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            one = root / "GCF" / "001" / "GCF_000001.1_genomic.fna.gz"
            two = root / "GCF" / "002" / "GCF_000002.1_genomic.fna.gz"
            one.parent.mkdir(parents=True)
            two.parent.mkdir(parents=True)
            with gzip.open(one, "wt", encoding="utf-8") as handle:
                handle.write(">a\nACGT\n")
            with gzip.open(two, "wt", encoding="utf-8") as handle:
                handle.write(">b\nACGT\n")

            self.assertEqual(discover_genome_files(root), [one, two])

    def test_output_paths_for_genome_mirrors_input_tree(self) -> None:
        genome = Path("00_raw_gtdb_r232/genomes/GCF/033/239/625/GCF_033239625.1_genomic.fna.gz")

        faa_path, gff_path = output_paths_for_genome(genome, Path("00_raw_gtdb_r232/genomes"), Path("03_gtdb_proteomes"))

        self.assertEqual(faa_path.as_posix(), "03_gtdb_proteomes/faa/GCF/033/239/625/GCF_033239625.1.faa.gz")
        self.assertEqual(gff_path.as_posix(), "03_gtdb_proteomes/gff/GCF/033/239/625/GCF_033239625.1.gff.gz")

    def test_write_prediction_outputs_uses_reversible_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            faa_path = root / "faa" / "GCF_000001.1.faa.gz"
            gff_path = root / "gff" / "GCF_000001.1.gff.gz"
            genome_result = {
                "accession": "GCF_000001.1",
                "genome_path": "source/GCF_000001.1_genomic.fna.gz",
                "contigs": [
                    {
                        "contig_id": "contigA",
                        "sequence": "ATG" + "AAA" * 30 + "TAA",
                        "genes": [
                            {
                                "begin": 1,
                                "end": 96,
                                "strand": 1,
                                "orf_index": 1,
                                "translation": "M" + "K" * 30,
                                "partial_begin": False,
                                "partial_end": False,
                            }
                        ],
                    }
                ],
            }

            write_prediction_outputs(genome_result, faa_path, gff_path)

            with gzip.open(faa_path, "rt", encoding="utf-8") as handle:
                faa_text = handle.read()
            with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
                gff_text = handle.read()

            self.assertIn(">GCF_000001.1|contigA|1", faa_text)
            self.assertIn("GCF_000001.1|contigA|1", gff_text)
            self.assertIn("ID=GCF_000001.1|contigA|1", gff_text)

    def test_predict_genome_uses_predictor_and_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genome_path = root / "GCF_000001.1_genomic.fna.gz"
            with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
                handle.write(">contigA\n")
                handle.write("ATG" + "AAA" * 30 + "TAA\n")

            policy = {"selected_predictor": "pyrodigal", "selected_mode": "meta"}

            calls: list[tuple[Path, str]] = []

            def fake_predictor(path: Path, accession: str, mode: str) -> dict[str, object]:
                calls.append((path, mode))
                return {
                    "accession": accession,
                    "genome_path": str(path),
                    "contigs": [
                        {
                            "contig_id": "contigA",
                            "sequence": "ATG" + "AAA" * 30 + "TAA",
                            "genes": [
                                {
                                    "begin": 1,
                                    "end": 96,
                                    "strand": 1,
                                    "orf_index": 1,
                                    "translation": "M" + "K" * 30,
                                    "partial_begin": False,
                                    "partial_end": False,
                                }
                            ],
                        }
                    ],
                    "summary": {
                        "contig_count": 1,
                        "predicted_genes": 1,
                        "total_bases": 96,
                        "coding_bases": 96,
                        "coding_density": 1.0,
                        "mean_protein_length": 31.0,
                        "internal_stops": 0,
                        "illegal_amino_acids": 0,
                        "short_orfs": 0,
                        "overlaps": 0,
                    },
                }

            result = predict_genome(genome_path, policy, predictor=fake_predictor)

            self.assertEqual(calls, [(genome_path, "meta")])
            self.assertEqual(result["accession"], "GCF_000001.1")
            self.assertEqual(result["summary"]["predicted_genes"], 1)

    def test_run_p03_prediction_preserves_input_order_with_threads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            g1 = root / "GCF_000001.1_genomic.fna.gz"
            g2 = root / "GCF_000002.1_genomic.fna.gz"
            for path in [g1, g2]:
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.write(">contigA\nATG" + "AAA" * 30 + "TAA\n")

            policy = {"selected_predictor": "pyrodigal", "selected_mode": "meta"}
            calls: list[str] = []

            def fake_predictor(path: Path, accession: str, mode: str) -> dict[str, object]:
                calls.append(accession)
                return {
                    "accession": accession,
                    "genome_path": str(path),
                    "contigs": [],
                    "summary": {
                        "contig_count": 1,
                        "predicted_genes": 1,
                        "total_bases": 96,
                        "coding_bases": 96,
                        "coding_density": 1.0,
                        "mean_protein_length": 31.0,
                        "internal_stops": 0,
                        "illegal_amino_acids": 0,
                        "short_orfs": 0,
                        "overlaps": 0,
                    },
                }

            result = run_p03_prediction(
                [g1, g2],
                policy,
                root / "out",
                threads=2,
                predictor=fake_predictor,
            )

            self.assertEqual(calls, ["GCF_000001.1", "GCF_000002.1"])
            self.assertEqual([row["accession"] for row in result["genomes"]], ["GCF_000001.1", "GCF_000002.1"])

    def test_load_prediction_policy_accepts_selected_predictor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            policy_path = Path(tmp) / "prediction_policy.yaml"
            policy_path.write_text(
                textwrap.dedent(
                    """
                    prediction:
                      stage: P02
                      selected_predictor: pyrodigal
                      selected_mode: meta
                      selected_interface: GeneFinder(meta=True)
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            policy = load_prediction_policy(policy_path)

            self.assertEqual(policy["selected_predictor"], "pyrodigal")
            self.assertEqual(policy["selected_mode"], "meta")


if __name__ == "__main__":
    unittest.main()
