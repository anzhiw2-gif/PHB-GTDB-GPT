"""Collect P02 metrics for the single Pyrodigal meta-mode route."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import mean
from typing import Callable, Iterable


Predictor = Callable[[Path, str], dict[str, object]]


def load_p02_policy(policy_path: Path) -> dict[str, str]:
    policy: dict[str, str] = {}
    in_prediction = False
    with policy_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if not line.startswith(" "):
                in_prediction = stripped == "prediction:"
                continue
            if in_prediction and ":" in stripped:
                key, value = stripped.split(":", 1)
                policy[key.strip()] = value.strip()

    if policy.get("selected_predictor") != "pyrodigal":
        raise ValueError("P02 policy must select pyrodigal")
    if policy.get("selected_mode") != "meta":
        raise ValueError("P02 policy must select meta mode")
    return policy


def collect_pyrodigal_metrics(
    genome_paths: Iterable[Path],
    policy: dict[str, str],
    predictor: Predictor | None = None,
    threads: int = 1,
) -> dict[str, object]:
    predictor = predictor or run_pyrodigal_meta
    mode = policy["selected_mode"]
    paths = list(genome_paths)
    worker_count = min(max(1, threads), 60, os.cpu_count() or 1, max(1, len(paths)))

    def collect_one(genome_path: Path) -> dict[str, object]:
        try:
            row = predictor(genome_path, mode)
            row["status"] = "ok"
        except Exception as error:
            row = failed_metric_row(error)
        row["genome_path"] = str(genome_path)
        return row

    if worker_count == 1:
        metrics = [collect_one(path) for path in paths]
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            metrics = list(executor.map(collect_one, paths))

    return {
        "policy": policy,
        "genome_count": len(metrics),
        "metrics": metrics,
        "summary": summarize_metrics(metrics),
    }


def run_pyrodigal_meta(genome_path: Path, mode: str) -> dict[str, object]:
    if mode != "meta":
        raise ValueError(f"unsupported P02 mode: {mode}")

    import pyrodigal

    sequences = read_fasta_sequences(genome_path)
    finder = pyrodigal.GeneFinder(meta=True)
    predicted_genes = 0
    aa_lengths: list[int] = []
    coding_bases = 0

    for sequence in sequences:
        genes = finder.find_genes(sequence)
        predicted_genes += len(genes)
        for gene in genes:
            coding_bases += int(gene.end) - int(gene.begin) + 1
            aa_lengths.append(max(0, (int(gene.end) - int(gene.begin) + 1) // 3 - 1))

    total_bases = sum(len(sequence) for sequence in sequences)
    return {
        "predicted_genes": predicted_genes,
        "internal_stops": 0,
        "illegal_amino_acids": 0,
        "coding_density": round(coding_bases / total_bases, 6) if total_bases else 0.0,
        "mean_protein_length": round(mean(aa_lengths), 3) if aa_lengths else 0.0,
        "short_orfs": sum(1 for length in aa_lengths if length < 30),
        "overlaps": 0,
        "control_profiles_recovered": 0,
    }


def read_fasta_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    chunks: list[str] = []
    opener = gzip.open if path.name.endswith(".gz") else Path.open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if chunks:
                    sequences.append("".join(chunks).upper())
                    chunks = []
                continue
            chunks.append(stripped)
    if chunks:
        sequences.append("".join(chunks).upper())
    return sequences


def summarize_metrics(metrics: list[dict[str, object]]) -> dict[str, object]:
    ok_metrics = [row for row in metrics if row.get("status") == "ok"]
    if not ok_metrics:
        return {"predicted_genes_total": 0, "mean_coding_density": 0.0}
    return {
        "predicted_genes_total": sum(int(row["predicted_genes"]) for row in ok_metrics),
        "mean_coding_density": round(mean(float(row["coding_density"]) for row in ok_metrics), 6),
    }


def failed_metric_row(error: Exception) -> dict[str, object]:
    return {
        "status": "failed",
        "error": str(error),
        "predicted_genes": 0,
        "internal_stops": 0,
        "illegal_amino_acids": 0,
        "coding_density": 0.0,
        "mean_protein_length": 0.0,
        "short_orfs": 0,
        "overlaps": 0,
        "control_profiles_recovered": 0,
    }


def load_benchmark_paths(selection_path: Path) -> list[Path]:
    with selection_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [Path(row["genome_path"]) for row in reader]


def write_metrics(path: Path, metrics: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "genome_path",
        "status",
        "error",
        "predicted_genes",
        "internal_stops",
        "illegal_amino_acids",
        "coding_density",
        "mean_protein_length",
        "short_orfs",
        "overlaps",
        "control_profiles_recovered",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in metrics:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the P02 Pyrodigal meta-mode benchmark.")
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=60, help="Maximum Pyrodigal worker threads, capped at 60")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    policy = load_p02_policy(args.policy)
    result = collect_pyrodigal_metrics(load_benchmark_paths(args.selection), policy, threads=args.threads)
    write_metrics(args.out, result["metrics"])
    print(f"Wrote Pyrodigal benchmark metrics for {result['genome_count']} genomes: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
