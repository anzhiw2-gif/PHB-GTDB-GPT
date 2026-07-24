from __future__ import annotations

import tempfile
import textwrap
import gzip
from pathlib import Path
import unittest

from scripts.p02_compare_predictors import collect_pyrodigal_metrics, load_p02_policy, read_fasta_sequences
from scripts.p02_select_benchmark_genomes import infer_accession, load_taxonomy, select_benchmark_genomes


class P02PredictorMetricsTest(unittest.TestCase):
    def test_select_benchmark_genomes_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genomes_dir = root / "genomes"
            genomes_dir.mkdir()
            taxonomy_path = root / "taxonomy.tsv"
            taxonomy_path.write_text(
                textwrap.dedent(
                    """
                    RS_GCF_000001.1\td__Bacteria;p__Firmicutes;c__Bacilli
                    RS_GCF_000002.1\td__Bacteria;p__Firmicutes;c__Bacilli
                    RS_GCF_000003.1\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria
                    GB_GCA_000004.1\td__Archaea;p__Euryarchaeota;c__Methanobacteria
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            for accession in ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCA_000004.1"]:
                genome_dir = genomes_dir / accession[:3] / accession
                genome_dir.mkdir(parents=True)
                (genome_dir / f"{accession}_genomic.fna").write_text(">contig\nACGT\n", encoding="utf-8")

            selected = select_benchmark_genomes(
                genomes_root=genomes_dir,
                taxonomy_path=taxonomy_path,
                sample_size=3,
                seed=20260724,
            )

            self.assertEqual([row["accession"] for row in selected], ["GCA_000004.1", "GCF_000001.1", "GCF_000003.1"])

    def test_infer_accession_uses_real_gtdb_filename(self) -> None:
        path = Path("genomes/GCF/033/239/625/GCF_033239625.1_genomic.fna.gz")

        self.assertEqual(infer_accession(path, Path("genomes")), "GCF_033239625.1")

    def test_collect_pyrodigal_metrics_reports_expected_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta_path = root / "sample.fna"
            fasta_path.write_text(">contig\nACGTACGTACGT\n", encoding="utf-8")
            policy_path = root / "prediction_policy.yaml"
            policy_path.write_text(
                "prediction:\n  selected_predictor: pyrodigal\n  selected_mode: meta\n",
                encoding="utf-8",
            )

            calls: list[tuple[Path, str]] = []

            def fake_predictor(path: Path, mode: str) -> dict[str, object]:
                calls.append((path, mode))
                return {
                    "predicted_genes": 2,
                    "internal_stops": 0,
                    "illegal_amino_acids": 0,
                    "coding_density": 0.84,
                    "mean_protein_length": 120.5,
                    "short_orfs": 0,
                    "overlaps": 0,
                    "control_profiles_recovered": 1,
                }

            result = collect_pyrodigal_metrics(
                genome_paths=[fasta_path],
                policy=load_p02_policy(policy_path),
                predictor=fake_predictor,
            )

            self.assertEqual(calls, [(fasta_path, "meta")])
            self.assertEqual(result["genome_count"], 1)
            self.assertEqual(result["metrics"][0]["predicted_genes"], 2)
            self.assertEqual(result["policy"]["selected_predictor"], "pyrodigal")

    def test_read_fasta_sequences_handles_gzip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = Path(tmp) / "sample.fna.gz"
            with gzip.open(fasta_path, "wt", encoding="utf-8") as handle:
                handle.write(">a\nacgt\n>b\nTTAA\n")

            self.assertEqual(read_fasta_sequences(fasta_path), ["ACGT", "TTAA"])

    def test_load_taxonomy_accepts_multiple_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bac = root / "bac.tsv"
            ar = root / "ar.tsv"
            bac.write_text("RS_GCF_000001.1\td__Bacteria;p__Firmicutes\n", encoding="utf-8")
            ar.write_text("GB_GCA_000002.1\td__Archaea;p__Methanobacteriota\n", encoding="utf-8")

            taxonomy = load_taxonomy([bac, ar])

            self.assertEqual(taxonomy["GCF_000001.1"], "d__Bacteria;p__Firmicutes")
            self.assertEqual(taxonomy["GCA_000002.1"], "d__Archaea;p__Methanobacteriota")

    def test_collect_pyrodigal_metrics_records_failed_genomes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good = root / "good.fna"
            bad = root / "bad.fna"
            good.write_text(">good\nACGT\n", encoding="utf-8")
            bad.write_text(">bad\nACGT\n", encoding="utf-8")
            policy = {"selected_predictor": "pyrodigal", "selected_mode": "meta"}

            def fake_predictor(path: Path, mode: str) -> dict[str, object]:
                if path.name == "bad.fna":
                    raise RuntimeError("cannot parse")
                return {
                    "predicted_genes": 1,
                    "internal_stops": 0,
                    "illegal_amino_acids": 0,
                    "coding_density": 0.5,
                    "mean_protein_length": 10.0,
                    "short_orfs": 1,
                    "overlaps": 0,
                    "control_profiles_recovered": 0,
                }

            result = collect_pyrodigal_metrics([good, bad], policy, predictor=fake_predictor)

            self.assertEqual(result["genome_count"], 2)
            self.assertEqual(result["metrics"][0]["status"], "ok")
            self.assertEqual(result["metrics"][1]["status"], "failed")
            self.assertEqual(result["metrics"][1]["error"], "cannot parse")

    def test_collect_pyrodigal_metrics_preserves_order_with_threads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = [root / "a.fna", root / "b.fna"]
            for path in paths:
                path.write_text(">x\nACGT\n", encoding="utf-8")

            def fake_predictor(path: Path, mode: str) -> dict[str, object]:
                return {
                    "predicted_genes": 1 if path.name == "a.fna" else 2,
                    "internal_stops": 0,
                    "illegal_amino_acids": 0,
                    "coding_density": 0.5,
                    "mean_protein_length": 10.0,
                    "short_orfs": 1,
                    "overlaps": 0,
                    "control_profiles_recovered": 0,
                }

            result = collect_pyrodigal_metrics(
                paths,
                {"selected_predictor": "pyrodigal", "selected_mode": "meta"},
                predictor=fake_predictor,
                threads=2,
            )

            self.assertEqual([row["genome_path"] for row in result["metrics"]], [str(paths[0]), str(paths[1])])


if __name__ == "__main__":
    unittest.main()
